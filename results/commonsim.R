
source("../src/simsel.R")
source("./defaults.R")
source("./studycases.R")
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
	theta           = rep(NA, default.n),
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

sim.run.single <- function(W0, myargs=NULL, sim.name=NA, force.run=FALSE) {
	missing.args <- names(default.args)[!(names(default.args) %in% names(myargs))]
	myargs[missing.args] <- default.args[missing.args]
	# Here the algorithm gets a bit complex, the problem is due to the fact that the number of genes can be deduced from different sources.
	# If W0 is NA, the W matrix needs to be built from another source of information (the length of theta). 
	# If W0 is a matrix, then we assume that the number of columns give the number of genes, and theta is updated accordingly
	# If both W0 and theta are NA, then the network is assumed to have only one gene. 
	if (length(W0) > 1) {
		if (length(myargs$theta) > ncol(W0)) myargs$theta <- myargs$theta[1:ncol(W0)]
		if (length(myargs$theta) < ncol(W0)) myargs$theta <- c(myargs$theta, rep(NA, ncol(W0)-length(myargs$theta)))
	}
	# NAs in the optimum are replaced by random values
	myargs$theta[is.na(myargs$theta)] <- runif(sum(is.na(myargs$theta)))
	if (all(is.na(W0))) {
		# NA for W0 is replaced by a random matrix at the correct equilibrium
		WW <- matrix(rnorm(length(myargs$theta)^2, default.rand.mean, default.rand.sd.sim), ncol=length(myargs$theta))
		for (i in 1:nrow(WW)) WW[i,sample(1:ncol(WW), 1)] <- NA
		W0 <- targetW(WW, target=myargs$theta, a=myargs$a)
	}
	
	if (is.na(sim.name)) 
		sim.name <- paste0("pure", paste(sample(c(letters, LETTERS), 10, replace=TRUE), collapse=""))
	outfile <- paste0(cache.dir, "/", sim.name ,".rds")
	if (force.run || !file.exists(outfile)) {
		ans <- do.call(simsel, c(list(W0=W0), myargs))
		saveRDS(ans, file=outfile, version=2)
	} else {
		ans <- readRDS(outfile)
	}
	ans
}

sim.run.reps <- function(W0, args=NULL, reps=10, series.name="pure", force.run=FALSE, mc.cores=detectCores()-1) {
	ans <- mclapply(1:reps, function(r) {
		ss <- sim.run.single(W0=W0, args, paste0(series.name, "-", r), force.run=force.run)
	}, mc.cores=mc.cores)

	list(full=ans, mean=acrossrepMean(ans), var=acrossrepVar(ans))
}


