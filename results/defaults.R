# 'source' this file before running any figure script

# Computing options

library(parallel)

default <- list(
	mc.cores = min(detectCores()-1, 100),

	use.cache = TRUE,
	cache.dir = "../cache/",
	pointsize = 6,   # For some reason, fonts appear larger than this setting
	maxfigwidth= 16, # 16 cm wide figures max
	figscale   =2.54, # cm to inches

	# General parameters

	mut.correlated = TRUE,
	a              = 0.2,
	n              = 6,
	dev.steps      = 16,
	dev.measure    = 4,
	epsilon.zero   = 0.01,

	# Random networks

	density        = 1,
	rand.mean      = 0,
	rand.sd        = 1,
	rand.sd.sim    = 1e-4,

	# Robustness index parameters

	rob.reps      = 100,
	initenv.sd    = 0.1,
	lateenv.sd    = 0.1,
	initmut.sd    = 0.1,
	latemut.sd    = 0.1,
	log.robustness= TRUE,
	summary.FUN   = mean,

	# Simulation parameters

	sim.reps      = 100,
	G             = 100,
	summary.every = 10,
	N             = 1000,
	initsd        = 0.001,
	mut.rate      = 0.01,
	sim.mutsd     = 0.1,
	s             = 10,
	nsel          = 3
)

# This has to be called somewhere before the simulation scripts, perhaps not the best here. 
if (!dir.exists(default$cache.dir)) dir.create(default$cache.dir)
