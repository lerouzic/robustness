
source("../src/puresel.R")
source("./defaults.R")
cache.dir <- "../cache"

library(parallel)

default.args <- list(
	a               = default.a,
	mut.correlated  = default.mut.correlated,
	dev.steps       = default.dev.steps,
	measure         = default.dev.measure,
	G               = default.G,
	summary.every   = round(default.G/100),
	mc.cores        = 1,
	N               = default.N,
	mut.rate        = default.mut.rate,
	som.mut.rate    = 0,
	sim.initenv.sd  = 0,
	sim.lateenv.sd  = 0,
	mut.sd          = default.sim.mutsd,
	theta           = 0.5,
	s               = default.s,
	grad.rob        = rep(0, 5),
	rep             = default.rob.reps,
	initenv.sd      = default.initenv.sd,
	lateenv.sd      = default.lateenv.sd,
	initmut.sd      = default.initmut.sd,
	latemut.sd      = default.latemut.sd,
	log.robustness  = default.log.robustness,
	plasticity      = FALSE
)
	
acrossrepMean <- function(fulllist) {
	ans <- lapply(names(fulllist[[1]]), function(gen) meanlist(lapply(fulllist, "[[", gen)))
	names(ans) <- names(fulllist[[1]])
	ans
}

acrossrepVar <- function(fulllist) {
	ans <- lapply(names(fulllist[[1]]), function(gen) varlist(lapply(fulllist, "[[", gen)))
	names(ans) <- names(fulllist[[1]])
	ans
}

pure.run.single <- function(W0, args=default.args, sim.name=NA, force.run=FALSE) {
	myargs <- default.args
	myargs[names(args)] <- args # Mistakes in arg names etc. will make the simulation crash later
	if (is.na(sim.name)) 
		sim.name <- paste0("pure", paste(sample(c(letters, LETTERS), 10, replace=TRUE), collapse=""))
	outfile <- paste0(cache.dir, "/", sim.name ,".rds")
	if (force.run || !file.exists(outfile)) {
		ans <- do.call(puresel, c(list(W0=W0), myargs))
		saveRDS(ans, file=outfile, version=2)
	} else {
		ans <- readRDS(outfile)
	}
	ans
}

pure.run.reps <- function(W0, args=NULL, reps=10, series.name="pure", force.run=FALSE, mc.cores=detectCores()-1) {
	ans <- mclapply(1:reps, function(r) {
		ss <- pure.run.single(W0=W0, args, paste0(series.name, "-", r), force.run=force.run)
	}, mc.cores=mc.cores)

	list(full=ans, mean=acrossrepMean(ans), var=acrossrepVar(ans))
}


