# =============================================================================
# example.R -- GCF demo on the paper's simulation data (900 grid cells)
#
# Pipeline: read data -> gcf_field() (generate GCF variables) ->
# gcf_select() (select a stable subset) -> print the selected variables.
# All settings below are the paper's simulation settings; the run takes
# about a minute.
# =============================================================================

source("gcf.R")

# Simulation data: response y1, covariates x1-x3, grid coordinates x, y
# (30 x 30 grid, 900 locations)
sim <- read.csv("sim-data.csv", row.names = 1)

# Step 1-3a: generate the GCF candidate variables (no response involved).
# Buffers 2, 4, 6 grid units; 11 quantile levels; LISA normalization radius
# d_norm = 4; fine/broad scale bands {2}/{6}.
field <- gcf_field(sim, coords = c("x", "y"), vars = c("x1", "x2", "x3"),
                   buffers = c(2, 4, 6), probs = seq(0, 1, 0.1),
                   d_norm = 4, fine_band = 2, broad_band = 6)
print(field)

# Step 3b: select a stable subset with spatial-block stability resampling.
# Blocks of 6 x 6 grid units; selection defaults are the paper settings
# (B = 80 resamples, ktop = 20, pi_thr = 0.6, seed = 1).
blocks <- gcf_blocks(sim[, c("x", "y")], size = 6)
sel <- gcf_select(field, y = sim$y1, blocks = blocks)
print(sel)

# The selected variables (raw covariates plus stable derived variables) are
# ready to feed any downstream regression learner, e.g. a random forest.
cat("\nSelected variables:\n")
print(sel$selected)
