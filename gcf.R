# =============================================================================
# gcf.R -- Generalized Covariate Field (GCF) model
#
# Author: Yongze Song
# Date:   29 August 2026
#
# Standalone implementation of the generalized covariate field (GCF) method:
# GCF variable generation (spatial-pattern features psi, neighbourhood-
# distribution features Zx, functional reduction) and variable selection
# (random forest importance + spatial-block stability + group voting).
# Consolidated from the R package `gcf` 0.1.0; function names and signatures
# are identical to the package.
#
# Dependencies (install once): sf, spdep, geocomplexity, ranger
# Requires R >= 4.1.
#
# Usage: source("gcf.R"), then see example.R. Exported functions:
#   gcf_field()  -- main generator: spatial variables in, GCF variables out
#   gcf_psi()    -- Step 1: spatial-pattern features (11 operators)
#   gcf_zx()     -- Step 2: neighbourhood-distribution (quantile) features
#   gcf_reduce() -- Step 3a: functional reduction to the candidate field
#   gcf_select() -- Step 3b: stable variable selection
#   gcf_blocks() -- square spatial block ids for stability/CV
# =============================================================================

library(sf)            # point geometry for the geocomplexity operator
library(spdep)         # k-nearest-neighbour weights for geocomplexity
library(geocomplexity) # geocomplexity operator (psi_G)
library(ranger)        # random forest importance kernel of gcf_select()

# =============================================================================
# Internal utilities: error helpers, input coercion, stage-1 cleaning
# =============================================================================

gcf_stop <- function(...) {
  stop(paste0("[gcf] ", paste0(..., collapse = "")), call. = FALSE)
}

gcf_assert <- function(ok, ...) {
  if (!isTRUE(ok)) gcf_stop(...)
  invisible(TRUE)
}

# Convert user data to a numeric data frame with stable names.
gcf_as_numeric_df <- function(X, arg_name = "X") {
  if (is.null(X)) return(NULL)
  X <- as.data.frame(X, stringsAsFactors = FALSE)
  if (!ncol(X)) gcf_stop(arg_name, " must have at least one column.")
  if (is.null(names(X)) || any(!nzchar(names(X)))) {
    gcf_stop(arg_name, " must have non-empty column names.")
  }
  if (anyDuplicated(names(X))) {
    gcf_stop(arg_name, " has duplicated column names: ",
             paste(unique(names(X)[duplicated(names(X))]), collapse = ", "))
  }
  for (nm in names(X)) {
    if (!is.numeric(X[[nm]])) {
      X[[nm]] <- suppressWarnings(as.numeric(X[[nm]]))
    }
    if (!is.numeric(X[[nm]])) gcf_stop(arg_name, "$", nm, " is not numeric.")
  }
  X
}

# Convert coordinates to a two-column numeric data frame.
gcf_as_coords <- function(coords) {
  coords <- as.data.frame(coords, stringsAsFactors = FALSE)
  gcf_assert(ncol(coords) >= 2, "coords must have at least two columns.")
  coords <- coords[, seq_len(2), drop = FALSE]
  names(coords) <- c("x", "y")
  coords$x <- suppressWarnings(as.numeric(coords$x))
  coords$y <- suppressWarnings(as.numeric(coords$y))
  gcf_assert(all(is.finite(coords$x)) && all(is.finite(coords$y)),
             "coords must be finite numeric values.")
  coords
}

# Stage-1 raw covariate cleaning: impute non-finite values by the column
# median (0 if the median is not finite) and drop zero-variance columns.
# Pattern and neighbourhood-distribution features are always built from this
# cleaned, raw-scale matrix.
gcf_raw_clean <- function(X) {
  X <- gcf_as_numeric_df(X, "X")
  impute_values <- vapply(X, function(col) {
    ok <- is.finite(col)
    med <- stats::median(col[ok], na.rm = TRUE)
    if (is.finite(med)) med else 0
  }, numeric(1))
  X_imp <- X
  for (nm in names(X_imp)) {
    bad <- !is.finite(X_imp[[nm]])
    if (any(bad)) X_imp[[nm]][bad] <- impute_values[[nm]]
  }
  zero_var <- vapply(X_imp, function(col) {
    vals <- unique(col[is.finite(col)])
    length(vals) <= 1L
  }, logical(1))
  retained <- names(X_imp)[!zero_var]
  dropped <- names(X_imp)[zero_var]
  gcf_assert(length(retained) > 0L,
             "Stage-1 cleaning dropped all covariate columns.")
  list(
    X = X_imp[, retained, drop = FALSE],
    retained = retained,
    dropped = dropped,
    impute_values = impute_values
  )
}

# Resolve the (data, coords, vars) user contract shared by gcf_psi, gcf_zx,
# and gcf_field: `coords` is either a length-2 character vector naming the
# coordinate columns of `data`, or a two-column matrix / data frame of
# projected coordinates; `vars` defaults to every other column of `data`.
gcf_resolve_input <- function(data, coords, vars = NULL) {
  gcf_assert(is.data.frame(data) || is.matrix(data),
             "data must be a data frame (or matrix).")
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  coord_cols <- character(0)
  if (is.character(coords)) {
    gcf_assert(length(coords) == 2L,
               "coords must name exactly two coordinate columns.")
    missing <- setdiff(coords, names(data))
    gcf_assert(!length(missing), "data is missing coordinate column(s): ",
               paste(missing, collapse = ", "))
    coord_cols <- coords
    coords_df <- gcf_as_coords(data[, coords, drop = FALSE])
  } else {
    coords_df <- gcf_as_coords(coords)
    gcf_assert(nrow(coords_df) == nrow(data),
               "coords must have one row per row of data.")
  }
  if (is.null(vars)) {
    vars <- setdiff(names(data), coord_cols)
  } else {
    vars <- as.character(vars)
    missing <- setdiff(vars, names(data))
    gcf_assert(!length(missing), "data is missing variable column(s): ",
               paste(missing, collapse = ", "))
  }
  gcf_assert(length(vars) > 0L, "At least one spatial variable is required.")
  list(data = data, coords = coords_df, vars = vars)
}

# Pairwise Euclidean distance matrix.
gcf_pairwise_dist <- function(query_coords, support_coords = query_coords) {
  query_coords <- gcf_as_coords(query_coords)
  support_coords <- gcf_as_coords(support_coords)
  same <- nrow(query_coords) == nrow(support_coords) &&
    isTRUE(all.equal(query_coords, support_coords, check.attributes = FALSE))
  if (same) return(as.matrix(stats::dist(as.matrix(query_coords))))
  dx <- outer(query_coords$x, support_coords$x, "-")
  dy <- outer(query_coords$y, support_coords$y, "-")
  sqrt(dx * dx + dy * dy)
}

# =============================================================================
# Step 1 -- spatial-pattern features (psi, the P layer). Per covariate the 11
# spatial operators of the GCF method (LISA, local Geary's c, log local
# variance, rank-quantile entropy, geocomplexity, log scale-variance, local
# variogram exponent, and signed z- and MAD-outlier strengths) over the
# configured buffer radii. Never uses a response.
# =============================================================================

psi_fmt_scale <- function(x) {
  format(x, trim = TRUE, scientific = FALSE)
}

# Feature-name / group metadata for the psi layer.
psi_multiscale_group_map <- function(vars, buffers, d_norm, gc_k = 23L,
                                     include_vario_exp = FALSE,
                                     vario_buffers = numeric(0)) {
  if (!isTRUE(include_vario_exp)) vario_buffers <- numeric(0)
  rows <- list()
  add_row <- function(v, tag, suffix, buffer = NA_real_) {
    rows[[length(rows) + 1L]] <<- data.frame(
      feature_name = paste0(v, "_", suffix),
      base_variable = v,
      group_id = paste0("P_", v),
      feature_category = tag,
      feature_type = "P",
      buffer = buffer,
      prob = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  for (v in vars) {
    for (b in buffers) {
      btxt <- psi_fmt_scale(b)
      add_row(v, "lisa", paste0("lisa_b", btxt, "_n", psi_fmt_scale(d_norm)), b)
    }
    for (b in buffers) {
      add_row(v, "geary", paste0("geary_b", psi_fmt_scale(b)), b)
    }
    for (b in buffers) {
      add_row(v, "lvar", paste0("lvar_b", psi_fmt_scale(b)), b)
    }
    for (b in buffers) {
      add_row(v, "qentropy", paste0("qentropy_b", psi_fmt_scale(b)), b)
    }
    add_row(v, "gc", paste0("gc_k", gc_k))
    add_row(v, "scalevar",
            paste0("scalevar_buffers", psi_fmt_scale(min(buffers)),
                   "_", psi_fmt_scale(max(buffers))))
    for (b in vario_buffers) {
      add_row(v, "vario_exp", paste0("vario_exp_b", psi_fmt_scale(b)), b)
    }
    for (b in buffers) {
      add_row(v, "poutlier_z", paste0("poutlier_z_b", psi_fmt_scale(b)), b)
    }
    for (b in buffers) {
      add_row(v, "noutlier_z", paste0("noutlier_z_b", psi_fmt_scale(b)), b)
    }
    for (b in buffers) {
      add_row(v, "poutlier_mad", paste0("poutlier_mad_b", psi_fmt_scale(b)), b)
    }
    for (b in buffers) {
      add_row(v, "noutlier_mad", paste0("noutlier_mad_b", psi_fmt_scale(b)), b)
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Build neighbour indices from a distance matrix.
gcf_neighbors_from_dist <- function(Dmat, radius, include_self = FALSE) {
  lapply(seq_len(nrow(Dmat)), function(i) {
    idx <- which(Dmat[i, ] <= radius & (include_self | Dmat[i, ] > 0))
    as.integer(idx)
  })
}

# Local normalization used internally for LISA; not itself a pattern feature.
psi_local_norm <- function(x, nb_norm) {
  z <- numeric(length(x))
  for (i in seq_along(x)) {
    idx <- nb_norm[[i]]
    if (length(idx) >= 2L) {
      m <- mean(x[idx], na.rm = TRUE)
      s <- stats::sd(x[idx], na.rm = TRUE)
      z[i] <- if (is.finite(s) && s > 0) (x[i] - m) / s else 0
    }
  }
  z
}

# LISA on locally normalized values.
psi_lisa <- function(z, nb_local) {
  vapply(seq_along(z), function(i) {
    idx <- nb_local[[i]]
    if (!length(idx)) return(0)
    zi <- mean(z[idx], na.rm = TRUE)
    out <- z[i] * zi
    if (is.finite(out)) out else 0
  }, numeric(1))
}

# Local Geary's c.
psi_geary <- function(x, nb_local, normalize = TRUE) {
  out <- vapply(seq_along(x), function(i) {
    idx <- nb_local[[i]]
    if (!length(idx)) return(0)
    mean((x[i] - x[idx])^2, na.rm = TRUE)
  }, numeric(1))
  if (isTRUE(normalize)) {
    v <- stats::var(x, na.rm = TRUE)
    if (is.finite(v) && v > 0) out <- out / v
  }
  out
}

# Log local variance.
psi_lvar <- function(x, nb_local) {
  log1p(vapply(seq_along(x), function(i) {
    idx <- nb_local[[i]]
    if (length(idx) < 2L) return(0)
    v <- stats::var(x[idx], na.rm = TRUE)
    if (is.finite(v)) v else 0
  }, numeric(1)))
}

# Rank-binned quantile entropy.
psi_qentropy <- function(x, nb_local, bins = 10) {
  pct <- rank(x, ties.method = "average", na.last = "keep") / sum(!is.na(x))
  xb <- as.integer(cut(pct, breaks = seq(0, 1, length.out = bins + 1L),
                       include.lowest = TRUE, labels = FALSE))
  out <- vapply(seq_along(x), function(i) {
    idx <- nb_local[[i]]
    if (length(idx) < 2L) return(0)
    counts <- tabulate(xb[idx], nbins = bins)
    freq <- counts[counts > 0L] / length(idx)
    val <- -sum(freq * log(freq)) / log(bins)
    if (is.finite(val)) val else 0
  }, numeric(1))
  out
}

# Log scale-variance across grid aggregations.
psi_scalevar <- function(x, coords, scales) {
  coords <- gcf_as_coords(coords)
  sv_mat <- matrix(NA_real_, nrow = length(x), ncol = length(scales))
  for (j in seq_along(scales)) {
    s <- scales[[j]]
    key <- paste(floor(coords$x / s), floor(coords$y / s), sep = "_")
    gm <- tapply(x, key, mean, na.rm = TRUE)
    sv_mat[, j] <- gm[key]
  }
  vals <- apply(sv_mat, 1, function(row) {
    v <- stats::var(row, na.rm = TRUE)
    if (is.finite(v)) v else 0
  })
  log1p(vals)
}

# Positive/negative z-score outlier strengths.
psi_outlier_z <- function(x, nb_local, theta = 2) {
  pos <- neg <- numeric(length(x))
  for (i in seq_along(x)) {
    idx <- nb_local[[i]]
    if (length(idx) >= 2L) {
      xi <- x[idx]
      s <- stats::sd(xi, na.rm = TRUE)
      if (is.finite(s) && s > 0) {
        z <- (xi - mean(xi, na.rm = TRUE)) / s
        pos[i] <- sum(abs(z[z > theta]), na.rm = TRUE)
        neg[i] <- sum(abs(z[z < -theta]), na.rm = TRUE)
      }
    }
  }
  data.frame(poutlier_z = pos, noutlier_z = neg)
}

# Positive/negative MAD outlier strengths.
psi_outlier_mad <- function(x, nb_local, theta = 2) {
  pos <- neg <- numeric(length(x))
  mad_floor <- 1e-3 * stats::mad(x, na.rm = TRUE) + 1e-9
  for (i in seq_along(x)) {
    idx <- nb_local[[i]]
    if (length(idx) >= 2L) {
      xi <- x[idx]
      ctr <- stats::median(xi, na.rm = TRUE)
      scl <- max(stats::mad(xi, na.rm = TRUE), mad_floor)
      z <- (xi - ctr) / scl
      pos[i] <- sum(abs(z[z > theta]), na.rm = TRUE)
      neg[i] <- sum(abs(z[z < -theta]), na.rm = TRUE)
    }
  }
  data.frame(poutlier_mad = pos, noutlier_mad = neg)
}

# Coordinate-only geocomplexity plan (k-nearest-neighbour weights).
psi_geocomplexity_plan <- function(coords, k = 23) {
  coords <- gcf_as_coords(coords)
  k <- min(as.integer(k), nrow(coords) - 1L)
  if (k < 1L) return(list(k = 0L))
  sf_base <- sf::st_as_sf(coords, coords = c("x", "y"), crs = sf::NA_crs_)
  nb <- spdep::knn2nb(spdep::knearneigh(sf_base, k = k))
  weights <- spdep::nb2mat(nb, style = "W")
  list(geometry = sf::st_geometry(sf_base), crs = sf::st_crs(sf_base),
       weights = weights, k = k)
}

# Geocomplexity vector from a saved plan.
psi_geocomplexity <- function(plan, x, variable) {
  if (isTRUE(plan$k < 1L)) return(numeric(length(x)))
  sf_data <- sf::st_sf(stats::setNames(data.frame(x), variable),
                       geometry = plan$geometry, crs = plan$crs)
  gc <- geocomplexity::geocd_vector(sf_data, wt = plan$weights,
                                    method = "moran", normalize = TRUE)
  val <- gc[[paste0("GC_", variable)]]
  val[!is.finite(val)] <- 0
  val
}

# Index pairs matching the lower-triangle order of stats::dist.
psi_dist_pairs <- function(idx) {
  m <- length(idx)
  np <- m * (m - 1L) / 2L
  a <- integer(np)
  b <- integer(np)
  k <- 1L
  for (col in seq_len(m - 1L)) {
    for (row in seq.int(col + 1L, m)) {
      a[k] <- idx[row]
      b[k] <- idx[col]
      k <- k + 1L
    }
  }
  list(a = a, b = b)
}

# Deterministic local-variogram plan; dependency-light backend, no external
# variogram package involved.
psi_vario_plan <- function(coords, nb_local, radius, nlags = 4L) {
  coords <- gcf_as_coords(coords)
  edges <- seq(0, radius, length.out = nlags + 1L)
  lapply(seq_len(nrow(coords)), function(i) {
    idx <- c(i, nb_local[[i]])
    m <- length(idx)
    if (m < 6L) return(NULL)
    d <- as.numeric(stats::dist(as.matrix(coords[idx, , drop = FALSE])))
    ok <- d > 0
    d <- d[ok]
    if (length(d) < nlags) return(NULL)
    pairs <- psi_dist_pairs(idx)
    pairs$a <- pairs$a[ok]
    pairs$b <- pairs$b[ok]
    lab <- cut(d, breaks = edges, include.lowest = TRUE, labels = FALSE)
    keep <- !is.na(lab)
    if (sum(keep) < nlags) return(NULL)
    d <- d[keep]
    lab <- lab[keep]
    pairs$a <- pairs$a[keep]
    pairs$b <- pairs$b[keep]
    hh <- tapply(d, lab, mean)
    wt <- tapply(d, lab, length)
    list(a = pairs$a, b = pairs$b, lab = lab, hh = hh, wt = wt)
  })
}

# Local variogram exponent from a saved plan.
psi_vario_exp <- function(x, plan) {
  ex <- numeric(length(x))
  for (i in seq_along(plan)) {
    p <- plan[[i]]
    if (is.null(p)) next
    sv <- (x[p$a] - x[p$b])^2 / 2
    gh <- tapply(sv, p$lab, mean)
    hh <- p$hh[names(gh)]
    wt <- p$wt[names(gh)]
    keep <- is.finite(gh) & gh > 0 & is.finite(hh) & hh > 0
    if (sum(keep) >= 2L) {
      lg <- as.numeric(log(gh[keep] + 1e-9))
      lh <- as.numeric(log(hh[keep]))
      w <- as.numeric(wt[keep])
      W <- sum(w)
      mx <- sum(w * lh) / W
      my <- sum(w * lg) / W
      sxx <- sum(w * (lh - mx)^2)
      sxy <- sum(w * (lh - mx) * (lg - my))
      if (is.finite(sxx) && sxx > 0) ex[i] <- sxy / sxx
    }
  }
  ex
}

# Fit the psi plan: neighbour lists, geocomplexity plan, variogram plans, and
# feature/group metadata.
psi_fit_plan <- function(X_support, coords_support, psi_buffers, d_norm,
                         theta, bins, include_vario_exp, vario_buffers) {
  X_support <- gcf_as_numeric_df(X_support, "X_support")
  coords_support <- gcf_as_coords(coords_support)
  psi_buffers <- sort(unique(as.numeric(psi_buffers)))
  gcf_assert(length(psi_buffers) > 0L && all(is.finite(psi_buffers)) &&
               all(psi_buffers > 0),
             "buffers must be positive numeric values.")
  d_norm <- as.numeric(d_norm)
  gcf_assert(length(d_norm) == 1L && is.finite(d_norm) && d_norm > 0,
             "d_norm must be a positive number.")
  scales <- psi_buffers
  theta <- as.numeric(theta)
  bins <- as.integer(bins)
  include_vario_exp <- isTRUE(include_vario_exp)
  vario_buffers <- if (include_vario_exp) {
    if (is.null(vario_buffers)) psi_buffers else as.numeric(vario_buffers)
  } else {
    numeric(0)
  }
  vario_buffers <- sort(unique(vario_buffers))
  if (include_vario_exp) {
    gcf_assert(length(vario_buffers) > 0L && all(is.finite(vario_buffers)) &&
                 all(vario_buffers > 0),
               "vario_buffers must be positive numeric values when include_vario_exp is TRUE.")
  }
  Dmat <- gcf_pairwise_dist(coords_support)
  neighbor_buffers <- sort(unique(c(psi_buffers, vario_buffers)))
  buffer_keys <- paste0("b", psi_fmt_scale(neighbor_buffers))
  nb_by_buffer <- stats::setNames(lapply(neighbor_buffers, function(b) {
    gcf_neighbors_from_dist(Dmat, b, include_self = FALSE)
  }), buffer_keys)
  nb_norm <- gcf_neighbors_from_dist(Dmat, d_norm, include_self = FALSE)
  gc_plan <- psi_geocomplexity_plan(coords_support, k = 23L)
  vario_plans <- if (include_vario_exp) {
    vario_keys <- paste0("b", psi_fmt_scale(vario_buffers))
    stats::setNames(lapply(seq_along(vario_buffers), function(i) {
      psi_vario_plan(coords_support, nb_by_buffer[[vario_keys[[i]]]],
                     radius = vario_buffers[[i]], nlags = 4L)
    }), vario_keys)
  } else {
    list()
  }
  group_map <- psi_multiscale_group_map(names(X_support), psi_buffers, d_norm,
                                        gc_k = 23L,
                                        include_vario_exp = include_vario_exp,
                                        vario_buffers = vario_buffers)
  list(
    vars = names(X_support),
    support_coords = coords_support,
    d_norm = d_norm,
    psi_buffers = psi_buffers,
    include_vario_exp = include_vario_exp,
    vario_buffers = vario_buffers,
    scales = scales,
    theta = theta,
    bins = bins,
    nb_by_buffer = nb_by_buffer,
    nb_norm = nb_norm,
    geocomplexity_plan = gc_plan,
    vario_plans = vario_plans,
    feature_names = group_map$feature_name,
    group_map = group_map
  )
}

# Apply the psi plan to its support rows.
psi_apply_plan <- function(plan, X_support) {
  X_support <- gcf_as_numeric_df(X_support, "X_support")
  X_support <- X_support[, plan$vars, drop = FALSE]
  out <- matrix(NA_real_, nrow = nrow(X_support), ncol = length(plan$feature_names))
  colnames(out) <- plan$feature_names
  col <- 0L
  for (v in plan$vars) {
    x <- X_support[[v]]
    z <- psi_local_norm(x, plan$nb_norm)
    block <- list()
    for (b in plan$psi_buffers) {
      bkey <- paste0("b", psi_fmt_scale(b))
      btxt <- psi_fmt_scale(b)
      nb <- plan$nb_by_buffer[[bkey]]
      block[[paste0(v, "_lisa_b", btxt, "_n", psi_fmt_scale(plan$d_norm))]] <- psi_lisa(z, nb)
    }
    for (b in plan$psi_buffers) {
      bkey <- paste0("b", psi_fmt_scale(b))
      btxt <- psi_fmt_scale(b)
      block[[paste0(v, "_geary_b", btxt)]] <- psi_geary(x, plan$nb_by_buffer[[bkey]])
    }
    for (b in plan$psi_buffers) {
      bkey <- paste0("b", psi_fmt_scale(b))
      btxt <- psi_fmt_scale(b)
      block[[paste0(v, "_lvar_b", btxt)]] <- psi_lvar(x, plan$nb_by_buffer[[bkey]])
    }
    for (b in plan$psi_buffers) {
      bkey <- paste0("b", psi_fmt_scale(b))
      btxt <- psi_fmt_scale(b)
      block[[paste0(v, "_qentropy_b", btxt)]] <- psi_qentropy(x, plan$nb_by_buffer[[bkey]],
                                                              plan$bins)
    }
    block[[paste0(v, "_gc_k23")]] <- psi_geocomplexity(plan$geocomplexity_plan, x, v)
    block[[paste0(v, "_scalevar_buffers", psi_fmt_scale(min(plan$psi_buffers)),
                  "_", psi_fmt_scale(max(plan$psi_buffers)))]] <-
      psi_scalevar(x, plan$support_coords, plan$scales)
    if (isTRUE(plan$include_vario_exp)) {
      for (b in plan$vario_buffers) {
        bkey <- paste0("b", psi_fmt_scale(b))
        btxt <- psi_fmt_scale(b)
        block[[paste0(v, "_vario_exp_b", btxt)]] <-
          psi_vario_exp(x, plan$vario_plans[[bkey]])
      }
    }
    outlier_z_by_buffer <- lapply(plan$psi_buffers, function(b) {
      psi_outlier_z(x, plan$nb_by_buffer[[paste0("b", psi_fmt_scale(b))]], plan$theta)
    })
    outlier_mad_by_buffer <- lapply(plan$psi_buffers, function(b) {
      psi_outlier_mad(x, plan$nb_by_buffer[[paste0("b", psi_fmt_scale(b))]], plan$theta)
    })
    for (i in seq_along(plan$psi_buffers)) {
      btxt <- psi_fmt_scale(plan$psi_buffers[[i]])
      block[[paste0(v, "_poutlier_z_b", btxt)]] <- outlier_z_by_buffer[[i]]$poutlier_z
    }
    for (i in seq_along(plan$psi_buffers)) {
      btxt <- psi_fmt_scale(plan$psi_buffers[[i]])
      block[[paste0(v, "_noutlier_z_b", btxt)]] <- outlier_z_by_buffer[[i]]$noutlier_z
    }
    for (i in seq_along(plan$psi_buffers)) {
      btxt <- psi_fmt_scale(plan$psi_buffers[[i]])
      block[[paste0(v, "_poutlier_mad_b", btxt)]] <- outlier_mad_by_buffer[[i]]$poutlier_mad
    }
    for (i in seq_along(plan$psi_buffers)) {
      btxt <- psi_fmt_scale(plan$psi_buffers[[i]])
      block[[paste0(v, "_noutlier_mad_b", btxt)]] <- outlier_mad_by_buffer[[i]]$noutlier_mad
    }
    block <- as.data.frame(block, stringsAsFactors = FALSE, check.names = FALSE)
    rng <- seq.int(col + 1L, col + ncol(block))
    gcf_assert(identical(colnames(out)[rng], names(block)),
               "psi block column order mismatch for ", v)
    out[, rng] <- as.matrix(block)
    col <- col + ncol(block)
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

# gcf_psi(): Step 1 of the GCF method. For each input variable, computes the
# 11 spatial-pattern operators over a series of buffer radii. `coords` is a
# length-2 character vector of column names or a two-column matrix/data
# frame; `vars` defaults to all non-coordinate columns (exclude the response
# when relying on the default). Returns a "gcf_psi" object with $features,
# $map, $vars, $params.
gcf_psi <- function(data, coords, vars = NULL, buffers,
                    d_norm = max(buffers), theta = 2, bins = 10,
                    include_vario_exp = TRUE, vario_buffers = NULL) {
  inp <- gcf_resolve_input(data, coords, vars)
  clean <- gcf_raw_clean(inp$data[, inp$vars, drop = FALSE])
  plan <- psi_fit_plan(clean$X, inp$coords, psi_buffers = buffers,
                       d_norm = d_norm, theta = theta, bins = bins,
                       include_vario_exp = include_vario_exp,
                       vario_buffers = vario_buffers)
  features <- psi_apply_plan(plan, clean$X)
  out <- list(
    features = features,
    map = plan$group_map,
    vars = plan$vars,
    dropped_vars = clean$dropped,
    params = list(buffers = plan$psi_buffers, d_norm = plan$d_norm,
                  theta = plan$theta, bins = plan$bins,
                  include_vario_exp = plan$include_vario_exp,
                  vario_buffers = plan$vario_buffers, gc_k = 23L)
  )
  class(out) <- "gcf_psi"
  out
}

print.gcf_psi <- function(x, ...) {
  cat("GCF spatial-pattern features (psi)\n")
  cat("  locations: ", nrow(x$features), "\n", sep = "")
  cat("  variables: ", length(x$vars), " (",
      paste(utils::head(x$vars, 5), collapse = ", "),
      if (length(x$vars) > 5) ", ..." else "", ")\n", sep = "")
  cat("  operators: ", length(unique(x$map$feature_category)),
      " | buffers: ", paste(x$params$buffers, collapse = ", "), "\n", sep = "")
  cat("  features:  ", ncol(x$features), "\n", sep = "")
  invisible(x)
}

# =============================================================================
# Step 2 -- neighbourhood-distribution features (Zx, the D layer). For each
# location v, variable x and buffer b, summarise the distribution of x over
# the locations within b by quantile levels tau:
#   Z_x(v; b, tau) = Q_tau({ x(u) : u in N(v, b) }).
# The buffer neighbourhood includes the location itself. Never uses a
# response.
# =============================================================================

# Fit the Zx plan: buffers, quantile levels, feature names, and group map.
zx_fit_plan <- function(X_support, coords_support, buffers, probs) {
  X_support <- gcf_as_numeric_df(X_support, "X_support")
  coords_support <- gcf_as_coords(coords_support)
  buffers <- sort(unique(as.numeric(buffers)))
  probs <- sort(unique(as.numeric(probs)))
  gcf_assert(length(buffers) > 0L && all(is.finite(buffers)) && all(buffers > 0),
             "buffers must be positive numeric values.")
  gcf_assert(length(probs) > 0L && all(is.finite(probs)) && all(probs >= 0) && all(probs <= 1),
             "probs must be in [0, 1].")
  feature_names <- unlist(lapply(names(X_support), function(v) {
    unlist(lapply(buffers, function(b) {
      paste0(v, "_b", format(b, trim = TRUE, scientific = FALSE),
             "_q", format(probs, trim = TRUE, scientific = FALSE))
    }), use.names = FALSE)
  }), use.names = FALSE)
  group_map <- do.call(rbind, lapply(names(X_support), function(v) {
    do.call(rbind, lapply(buffers, function(b) {
      data.frame(buffer = b, prob = probs, stringsAsFactors = FALSE)
    })) |>
      transform(
        feature_name = paste0(v, "_b", format(buffer, trim = TRUE, scientific = FALSE),
                              "_q", format(prob, trim = TRUE, scientific = FALSE)),
        base_variable = v,
        group_id = paste0("D_", v),
        feature_category = "context_quantile",
        feature_type = "D"
      ) |>
      subset(select = c("feature_name", "base_variable", "group_id",
                        "feature_category", "feature_type", "buffer", "prob"))
  }))
  rownames(group_map) <- NULL
  list(
    vars = names(X_support),
    support_coords = coords_support,
    buffers = buffers,
    probs = probs,
    feature_names = feature_names,
    group_map = group_map
  )
}

# Apply the Zx plan. Empty buffers deliberately produce NA.
zx_apply_plan <- function(plan, X_support, coords_query,
                          query_to_support_row = NULL) {
  X_support <- gcf_as_numeric_df(X_support, "X_support")
  X_support <- X_support[, plan$vars, drop = FALSE]
  coords_query <- gcf_as_coords(coords_query)
  if (!is.null(query_to_support_row)) {
    gcf_assert(length(query_to_support_row) == nrow(coords_query),
               "query_to_support_row length must match coords_query rows.")
    Dmat <- gcf_pairwise_dist(plan$support_coords[query_to_support_row, , drop = FALSE],
                              plan$support_coords)
  } else {
    Dmat <- gcf_pairwise_dist(coords_query, plan$support_coords)
  }
  out <- matrix(NA_real_, nrow = nrow(coords_query), ncol = length(plan$feature_names))
  colnames(out) <- plan$feature_names
  col_offset <- 0L
  for (v in plan$vars) {
    x <- X_support[[v]]
    for (b in plan$buffers) {
      idx_by_row <- lapply(seq_len(nrow(Dmat)), function(i) which(Dmat[i, ] <= b))
      for (q in plan$probs) {
        col_offset <- col_offset + 1L
        out[, col_offset] <- vapply(idx_by_row, function(idx) {
          if (!length(idx)) return(NA_real_)
          as.numeric(stats::quantile(x[idx], probs = q, na.rm = TRUE, names = FALSE, type = 7))
        }, numeric(1))
      }
    }
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

# gcf_zx(): Step 2 of the GCF method. Buffer-wise quantiles of each variable
# over the neighbourhood of each location. `probs` defaults to
# seq(0, 1, 0.05) (21 levels, as in the paper's case study). Returns a
# "gcf_zx" object with $features (columns named <var>_b<buffer>_q<prob>),
# $map, $vars, $params.
gcf_zx <- function(data, coords, vars = NULL, buffers,
                   probs = seq(0, 1, 0.05)) {
  inp <- gcf_resolve_input(data, coords, vars)
  clean <- gcf_raw_clean(inp$data[, inp$vars, drop = FALSE])
  plan <- zx_fit_plan(clean$X, inp$coords, buffers = buffers, probs = probs)
  features <- zx_apply_plan(plan, clean$X, inp$coords,
                            query_to_support_row = seq_len(nrow(clean$X)))
  out <- list(
    features = features,
    map = plan$group_map,
    vars = plan$vars,
    dropped_vars = clean$dropped,
    params = list(buffers = plan$buffers, probs = plan$probs)
  )
  class(out) <- "gcf_zx"
  out
}

print.gcf_zx <- function(x, ...) {
  cat("GCF neighbourhood-distribution features (Zx)\n")
  cat("  locations: ", nrow(x$features), "\n", sep = "")
  cat("  variables: ", length(x$vars), " (",
      paste(utils::head(x$vars, 5), collapse = ", "),
      if (length(x$vars) > 5) ", ..." else "", ")\n", sep = "")
  cat("  buffers:   ", paste(x$params$buffers, collapse = ", "),
      " | quantile levels: ", length(x$params$probs), "\n", sep = "")
  cat("  features:  ", ncol(x$features), "\n", sep = "")
  invisible(x)
}

# =============================================================================
# Step 3a -- functional reduction. Collapses the collinear quantile/buffer
# sweeps of the psi and Zx layers into a compact set of interpretable
# functionals per (variable, scale band):
#   D layer: median, IQR, low tail (q0.10), high tail (q0.90), skew
#            (q0.90 + q0.10 - 2 q0.50) of the band-averaged quantile curve;
#   P layer: band-averaged buffered operators (single-scale gc and scalevar
#            kept as-is);
#   X layer: the raw covariates passed through unchanged.
# =============================================================================

# Core reduction: data holds raw covariate columns plus all psi and Zx
# feature columns; manifest has block ("X"/"P"/"D"), feature_name,
# base_variable, feature_category, buffer, prob.
gcf_build_reduced_candidates <- function(data, manifest,
                                         fine_band  = c(20, 30),
                                         broad_band = c(90, 100),
                                         vars = NULL,
                                         d_mode = c("functional",
                                                    "functional_allscale",
                                                    "qgrid", "full"),
                                         d_qgrid = c(0.05, 0.10, 0.25, 0.50,
                                                     0.75, 0.90, 0.95)) {
  d_mode <- match.arg(d_mode)
  man <- manifest
  if (is.null(vars)) vars <- unique(man$base_variable[man$block == "X"])

  bands <- list(fine = fine_band, broad = broad_band)
  cols  <- list()   # named numeric vectors
  meta  <- list()   # rows: feature, group, category, operator, scale

  add <- function(name, vec, group, category,
                  operator = NA_character_, scale = NA_character_) {
    cols[[name]]  <<- vec
    meta[[length(meta) + 1L]] <<- data.frame(
      feature = name, group = group, category = category,
      operator = operator, scale = scale, stringsAsFactors = FALSE)
  }

  # helper: band-averaged column for a manifest row subset
  band_avg <- function(colnames_in_band) {
    if (length(colnames_in_band) == 1L) return(data[[colnames_in_band]])
    rowMeans(as.matrix(data[, colnames_in_band, drop = FALSE]), na.rm = TRUE)
  }

  for (v in vars) {

    # ---- X: raw covariate (forced) ----
    add(v, data[[v]], group = v, category = "X")

    # ---- D: context-quantile features (representation set by d_mode) ----
    dman <- man[man$block == "D" & man$base_variable == v, ]
    if (d_mode %in% c("functional", "functional_allscale")) {
      # functionals of the quantile curve, per scale band
      d_bands <- if (d_mode == "functional") bands else
        stats::setNames(lapply(sort(unique(dman$buffer)), identity),
                        paste0("b", sort(unique(dman$buffer))))
      for (bn in names(d_bands)) {
        bset <- d_bands[[bn]]
        qcol <- function(p) {
          sub <- dman[dman$buffer %in% bset, ]
          avail <- unique(sub$prob)
          pn <- avail[which.min(abs(avail - p))]   # nearest available quantile
          band_avg(sub$feature_name[abs(sub$prob - pn) < 1e-9])
        }
        q10 <- qcol(0.10); q25 <- qcol(0.25); q50 <- qcol(0.50)
        q75 <- qcol(0.75); q90 <- qcol(0.90)
        pre <- paste0(v, "_D_", bn, "_")
        add(paste0(pre, "med"),    q50,                 v, "D", scale = bn)
        add(paste0(pre, "iqr"),    q75 - q25,           v, "D", scale = bn)
        add(paste0(pre, "lotail"), q10,                 v, "D", scale = bn)
        add(paste0(pre, "hitail"), q90,                 v, "D", scale = bn)
        add(paste0(pre, "skew"),   q90 + q10 - 2 * q50, v, "D", scale = bn)
      }
    } else if (d_mode == "qgrid") {
      # raw quantile values, all buffers x a coarse quantile grid (no functional)
      keep <- dman[vapply(dman$prob, function(p)
        any(abs(p - d_qgrid) < 1e-9), logical(1)), ]
      for (j in seq_len(nrow(keep))) {
        add(keep$feature_name[j], data[[keep$feature_name[j]]], v, "D",
            scale = paste0("b", keep$buffer[j]))
      }
    } else {  # "full": all D columns as-is, no reduction
      for (j in seq_len(nrow(dman))) {
        add(dman$feature_name[j], data[[dman$feature_name[j]]], v, "D",
            scale = paste0("b", dman$buffer[j]))
      }
    }

    # ---- P: pattern operators, band-averaged ----
    pman <- man[man$block == "P" & man$base_variable == v, ]
    ops  <- unique(pman$feature_category)
    for (op in ops) {
      opman <- pman[pman$feature_category == op, ]
      has_buffers <- any(!is.na(opman$buffer))
      if (has_buffers) {
        for (bn in names(bands)) {
          bset <- bands[[bn]]
          nm <- opman$feature_name[opman$buffer %in% bset]
          if (length(nm) == 0L) next
          add(paste0(v, "_P_", op, "_", bn), band_avg(nm), v, "P",
              operator = op, scale = bn)
        }
      } else {
        # single-scale operator (gc, scalevar)
        add(paste0(v, "_P_", op), data[[opman$feature_name[1]]], v, "P",
            operator = op, scale = "single")
      }
    }
  }

  meta_df <- do.call(rbind, meta)
  X <- do.call(cbind, cols)
  colnames(X) <- names(cols)
  rownames(meta_df) <- NULL
  stopifnot(identical(colnames(X), meta_df$feature))
  list(X = X, meta = meta_df)
}

# Assemble the manifest expected by gcf_build_reduced_candidates from the
# fitted psi/Zx layer metadata.
gcf_assemble_manifest <- function(psi_map, zx_map, vars) {
  keep <- c("block", "feature_name", "base_variable", "feature_category",
            "buffer", "prob")
  rbind(
    data.frame(block = "X", feature_name = vars, base_variable = vars,
               feature_category = "raw", buffer = NA_real_, prob = NA_real_,
               stringsAsFactors = FALSE),
    transform(psi_map, block = feature_type)[, keep],
    transform(zx_map, block = feature_type)[, keep]
  )
}

# gcf_reduce(): Step 3a of the GCF method. Reduces the psi and Zx sweeps to
# the compact candidate field over two scale bands (fine_band/broad_band
# default to min/max of the buffers used; the paper's case study uses
# c(20, 30) and c(90, 100) km). d_mode "functional" (default, as in the
# paper), "functional_allscale", "qgrid", or "full". Returns a "gcf_field"
# object with $candidates, $meta, $vars, $psi, $zx, $params.
gcf_reduce <- function(data, psi, zx, fine_band = NULL, broad_band = NULL,
                       d_mode = "functional") {
  gcf_assert(inherits(psi, "gcf_psi"), "psi must be a gcf_psi object.")
  gcf_assert(inherits(zx, "gcf_zx"), "zx must be a gcf_zx object.")
  gcf_assert(identical(psi$vars, zx$vars),
             "psi and zx must be built from the same variables.")
  gcf_assert(nrow(psi$features) == nrow(zx$features),
             "psi and zx must be built from the same locations.")
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  gcf_assert(nrow(data) == nrow(psi$features),
             "data must have the same rows as the psi/zx features.")
  missing <- setdiff(psi$vars, names(data))
  gcf_assert(!length(missing), "data is missing variable column(s): ",
             paste(missing, collapse = ", "))
  bufs <- sort(unique(c(psi$params$buffers, zx$params$buffers)))
  if (is.null(fine_band)) fine_band <- min(bufs)
  if (is.null(broad_band)) broad_band <- max(bufs)
  if (d_mode == "functional") {
    gcf_assert(any(zx$params$buffers %in% fine_band),
               "fine_band must contain at least one Zx buffer radius.")
    gcf_assert(any(zx$params$buffers %in% broad_band),
               "broad_band must contain at least one Zx buffer radius.")
  }
  modeling <- cbind(data[, psi$vars, drop = FALSE], psi$features, zx$features)
  man <- gcf_assemble_manifest(psi$map, zx$map, psi$vars)
  bad <- setdiff(man$feature_name, names(modeling))
  gcf_assert(!length(bad),
             "Feature name/manifest mismatch (buffer radii of mixed decimal ",
             "widths, e.g. c(1.5, 3), are not supported; use a consistent ",
             "buffer series): ", paste(utils::head(bad, 3), collapse = ", "))
  red <- gcf_build_reduced_candidates(modeling, man,
                                      fine_band = fine_band,
                                      broad_band = broad_band,
                                      d_mode = d_mode)
  out <- list(
    candidates = as.data.frame(red$X, stringsAsFactors = FALSE),
    meta = red$meta,
    vars = psi$vars,
    psi = psi,
    zx = zx,
    params = list(fine_band = fine_band, broad_band = broad_band,
                  d_mode = d_mode)
  )
  class(out) <- "gcf_field"
  out
}

# =============================================================================
# Main generator -- raw spatial variables in, their generalized covariate
# field variables out (psi -> Zx -> functional reduction).
# =============================================================================

# gcf_field(): the main GCF variable generator. Runs gcf_psi(), gcf_zx(), and
# gcf_reduce() in one call. No response variable is used at any point of the
# generation. Returns a "gcf_field" object: $candidates (raw X plus derived
# P and D variables), $meta (feature, group, category, operator, scale),
# $vars, $psi, $zx, $params.
gcf_field <- function(data, coords, vars = NULL, buffers,
                      probs = seq(0, 1, 0.05),
                      d_norm = max(buffers), theta = 2, bins = 10,
                      include_vario_exp = TRUE, vario_buffers = NULL,
                      fine_band = min(buffers), broad_band = max(buffers),
                      d_mode = "functional") {
  inp <- gcf_resolve_input(data, coords, vars)
  psi <- gcf_psi(inp$data, inp$coords, vars = inp$vars, buffers = buffers,
                 d_norm = d_norm, theta = theta, bins = bins,
                 include_vario_exp = include_vario_exp,
                 vario_buffers = vario_buffers)
  zx <- gcf_zx(inp$data, inp$coords, vars = inp$vars, buffers = buffers,
               probs = probs)
  gcf_reduce(inp$data, psi, zx, fine_band = fine_band,
             broad_band = broad_band, d_mode = d_mode)
}

print.gcf_field <- function(x, ...) {
  tab <- table(factor(x$meta$category, levels = c("X", "P", "D")))
  cat("Generalized covariate field (GCF)\n")
  cat("  locations:  ", nrow(x$candidates), "\n", sep = "")
  cat("  variables:  ", length(x$vars), " (",
      paste(utils::head(x$vars, 5), collapse = ", "),
      if (length(x$vars) > 5) ", ..." else "", ")\n", sep = "")
  cat("  candidates: ", ncol(x$candidates),
      "  [X (raw) ", tab[["X"]], " | P (pattern) ", tab[["P"]],
      " | D (context) ", tab[["D"]], "]\n", sep = "")
  invisible(x)
}

summary.gcf_field <- function(object, ...) {
  print(object)
  cat("\nCandidate variables per input variable and category:\n")
  tab <- table(object$meta$group, factor(object$meta$category,
                                         levels = c("X", "P", "D")))
  print(tab[object$vars, , drop = FALSE])
  cat("\nReduction: d_mode = ", object$params$d_mode,
      " | fine band {", paste(object$params$fine_band, collapse = ", "),
      "} | broad band {", paste(object$params$broad_band, collapse = ", "),
      "}\n", sep = "")
  cat("Full sweeps kept in $psi (", ncol(object$psi$features),
      " pattern features) and $zx (", ncol(object$zx$features),
      " context-quantile features).\n", sep = "")
  invisible(object)
}

# =============================================================================
# Step 3b -- variable selection. Random forest impurity importance (top-K per
# subsample) combined with spatial-block stability resampling and medium
# (variable x category) group voting.
# =============================================================================

# Core selection loop: X is the reduced candidate matrix, meta has
# feature/group/category, block_id covers all rows. Deterministic for a
# given seed (ranger runs single-threaded with a fixed seed).
gcf_select_rfimp_core <- function(X, y, meta, train_rows, block_id,
                                  params = list()) {
  p <- utils::modifyList(list(B = 80L, subsample_frac = 0.7, pi = 0.6,
                              ktop = 20L, rf_trees = 200L, seed = 1L), params)
  is_X <- meta$category == "X"; pen <- which(!is_X); x_cols <- meta$feature[is_X]
  set.seed(p$seed)
  Ztr <- X[train_rows, , drop = FALSE]; ytr <- y[train_rows]
  btr <- block_id[train_rows]; blocks <- unique(btr)
  kt <- floor(p$subsample_frac * length(blocks))
  sel <- matrix(FALSE, p$B, length(pen))
  for (b in seq_len(p$B)) {
    sub <- which(btr %in% sample(blocks, kt))
    rf <- ranger::ranger(x = Ztr[sub, , drop = FALSE], y = ytr[sub],
                         num.trees = p$rf_trees, importance = "impurity",
                         seed = 1, num.threads = 1)
    imp <- rf$variable.importance[meta$feature[pen]]; imp[is.na(imp)] <- 0
    sel[b, ] <- rank(-imp, ties.method = "first") <= p$ktop
  }
  key <- paste(meta$group[pen], meta$category[pen])
  reps <- unlist(lapply(unique(key), function(kk) {
    idx <- which(key == kk)
    if (mean(rowSums(sel[, idx, drop = FALSE]) > 0) < p$pi) return(NULL)
    meta$feature[pen][idx][which.max(colMeans(sel[, idx, drop = FALSE]))]
  }))
  list(selected = c(x_cols, reps), x_cols = x_cols, reps = reps,
       sel = sel, pen = pen, key = key, params = p)
}

# gcf_blocks(): assigns each location to a square spatial block of side
# `size` on projected coordinates. Block ids drive the stability resampling
# of gcf_select() and can also define spatial cross-validation folds.
gcf_blocks <- function(coords, size) {
  coords <- gcf_as_coords(coords)
  gcf_assert(length(size) == 1L && is.finite(size) && size > 0,
             "size must be a positive number.")
  cx <- coords[, 1]; cy <- coords[, 2]
  paste0("bx", floor((cx - min(cx)) / size), "_by", floor((cy - min(cy)) / size))
}

# gcf_select(): Step 3b of the GCF method. Screens the candidate variables of
# a gcf_field for a stable subset. The raw covariates (category "X") are
# always kept; the derived P/D variables are selected by (1) per-subsample
# random forest impurity importance (top `ktop`), (2) `B` spatial-block
# stability resamples drawing `subsample_frac` of the blocks, and (3) medium
# (variable x category) group voting with fire threshold `pi_thr` -- each
# qualifying group contributes its most frequently kept member. Defaults are
# the paper settings (B = 80, subsample_frac = 0.7, pi_thr = 0.6, ktop = 20,
# num_trees = 200, seed = 1). Returns a "gcf_selection" object with
# $selected, $forced, $derived, $freq, $group_fire, $params.
gcf_select <- function(x, y, blocks, train = NULL, B = 80,
                       subsample_frac = 0.7, pi_thr = 0.6, ktop = 20,
                       num_trees = 200, seed = 1) {
  gcf_assert(is.list(x) && !is.null(x$candidates) && !is.null(x$meta),
             "x must be a gcf_field object (or a list with candidates and meta).")
  X <- as.matrix(x$candidates)
  gcf_assert(is.numeric(X), "candidate variables must be numeric.")
  meta <- x$meta
  y <- as.numeric(y)
  gcf_assert(length(y) == nrow(X), "y must have one value per location.")
  gcf_assert(all(is.finite(y)), "y must be finite.")
  gcf_assert(length(blocks) == nrow(X),
             "blocks must have one id per location; see gcf_blocks().")
  if (is.null(train)) train <- seq_len(nrow(X))
  core <- gcf_select_rfimp_core(
    X, y, meta, train_rows = train, block_id = blocks,
    params = list(B = as.integer(B), subsample_frac = subsample_frac,
                  pi = pi_thr, ktop = as.integer(ktop),
                  rf_trees = as.integer(num_trees), seed = seed))
  freq <- colMeans(core$sel)
  names(freq) <- meta$feature[core$pen]
  group_fire <- vapply(unique(core$key), function(kk) {
    idx <- which(core$key == kk)
    mean(rowSums(core$sel[, idx, drop = FALSE]) > 0)
  }, numeric(1))
  out <- list(
    selected = core$selected,
    forced = core$x_cols,
    derived = if (is.null(core$reps)) character(0) else core$reps,
    freq = sort(freq, decreasing = TRUE),
    group_fire = sort(group_fire, decreasing = TRUE),
    params = list(B = as.integer(B), subsample_frac = subsample_frac,
                  pi_thr = pi_thr, ktop = as.integer(ktop),
                  num_trees = as.integer(num_trees), seed = seed,
                  n_train = length(train)),
    n_candidates = ncol(X)
  )
  class(out) <- "gcf_selection"
  out
}

print.gcf_selection <- function(x, ...) {
  cat("GCF variable selection (rf_imp + spatial-block stability + group voting)\n")
  cat("  candidates: ", x$n_candidates, " | resamples B = ", x$params$B,
      " | pi threshold = ", x$params$pi_thr, "\n", sep = "")
  cat("  selected:   ", length(x$selected), " (", length(x$forced),
      " forced raw + ", length(x$derived), " derived)\n", sep = "")
  cat("  forced:  ", paste(x$forced, collapse = ", "), "\n", sep = "")
  if (length(x$derived)) {
    cat("  derived: ", paste(x$derived, collapse = ", "), "\n", sep = "")
  } else {
    cat("  derived: (none passed the group vote)\n")
  }
  invisible(x)
}
