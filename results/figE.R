#!/usr/bin/env Rscript

source("./commonsim.R")

library(parallel)
mc.cores <- min(64, detectCores()-1)

# Standard stabilizing selection
ss <- sim.run(default.args, "try1", force.run=TRUE)
