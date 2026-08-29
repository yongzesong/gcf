# Generalized Covariate Field (GCF): Code and Examples

The generalized covariate field (GCF) model is a prediction-oriented method for spatial variables: it expands raw spatial covariates into spatial-pattern and neighbourhood-distribution features and selects a stable subset of them for geospatial prediction. For each covariate `x` observed at projected coordinates, the method first computes spatial-pattern features `psi` (11 spatial operators — LISA, local Geary's c, log local variance, rank-quantile entropy, geocomplexity, log scale-variance, local variogram exponent, and signed z-score and MAD outlier strengths — over a series of buffer radii), then neighbourhood-distribution features `Z_x` (buffer-wise quantiles of the covariate values surrounding each location), reduces the collinear buffer and quantile sweeps to a compact set of interpretable functionals over two scale bands, and finally selects variables by random forest importance combined with spatial-block stability resampling and group voting. No response variable is used at any point of the feature generation; the selected variables feed any downstream regression learner.

The pipeline: `x -> psi (pattern) -> Z_x (context) -> functional reduction -> selection`, implemented in `gcf.R` as:

1. `gcf_psi()` — Step 1: spatial-pattern features (11 operators over buffer radii);
2. `gcf_zx()` — Step 2: neighbourhood-distribution features (buffer-wise quantiles);
3. `gcf_reduce()` — Step 3a: functional reduction of the sweeps to the candidate field (X raw, P pattern, D context variables);
4. `gcf_select()` — Step 3b: stable variable selection (random forest importance + spatial-block stability + group voting), with the block helper `gcf_blocks()`;
5. `gcf_field()` — runs Steps 1–3a in one call.

## Usage

Dependencies (install once from CRAN): `sf`, `spdep`, `geocomplexity`, `ranger`. Requires R >= 4.1.

`example.R` runs the full pipeline on the paper's simulation data (`sim-data.csv`: 900 grid cells, response `y1`, covariates `x1`–`x3`, coordinates `x`, `y`) with the paper's settings, in about a minute:

``` r
source("gcf.R")

sim <- read.csv("sim-data.csv", row.names = 1)

# Generate the GCF candidate variables (no response involved)
field <- gcf_field(sim, coords = c("x", "y"), vars = c("x1", "x2", "x3"),
                   buffers = c(2, 4, 6), probs = seq(0, 1, 0.1),
                   d_norm = 4, fine_band = 2, broad_band = 6)

# Select a stable subset with spatial-block stability resampling
blocks <- gcf_blocks(sim[, c("x", "y")], size = 6)
sel <- gcf_select(field, y = sim$y1, blocks = blocks)
sel$selected
```

Output:

```
Generalized covariate field (GCF)
  locations:  900
  variables:  3 (x1, x2, x3)
  candidates: 93  [X (raw) 3 | P (pattern) 60 | D (context) 30]
GCF variable selection (rf_imp + spatial-block stability + group voting)
  candidates: 93 | resamples B = 80 | pi threshold = 0.6
  selected:   8 (3 forced raw + 5 derived)
  forced:  x1, x2, x3
  derived: x1_D_fine_med, x2_D_fine_med, x2_P_gc, x3_D_fine_med, x3_P_gc
```

The selection is deterministic for a given seed (`seed = 1` by default; the random forests run single-threaded with a fixed seed).
