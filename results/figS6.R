source("./terminology.R")
source("./defaults.R")
source("../src/robindex.R")

library(parallel)

a             <- default.a
dev.steps     <- default.dev.steps
dev.measure   <- default.dev.measure
log.robindex  <- default.log.robustness

rob.reps      <- 1000

net.size      <- default.n
reps <- 8

density        <- 1
reg.mean       <- -0.2
reg.sd         <- 1.2

mc.cores       <- default.mc.cores

nb.values <- 21
sd.test <- list(
	initenv = 10^seq(-4, 0, length.out=nb.values),
	lateenv = 10^seq(-4, 0, length.out=nb.values),
	initmut = 10^seq(-4, 0, length.out=nb.values),
	latemut = 10^seq(-4, 0, length.out=nb.values)
)

phen <- phen.expression
    
sigmas <- c(
	initenv=substitute(sigma[x], list(x=SDLETTER.ENVCAN)), 
	lateenv=substitute(sigma[x], list(x=SDLETTER.HOMEO)),
	initmut=substitute(sigma[x], list(x=SDLETTER.GENCAN)),
	latemut=substitute(sigma[x], list(x=SDLETTER.SOM))
)

col.phen <- c(initenv=COL.ENVCAN, lateenv=COL.HOMEO, initmut=COL.GENCAN, latemut=COL.SOM, stability=COL.STAB)

randW <- function(size, density, mean, sd) {
	W <- matrix(rnorm(size^2, mean=mean, sd=sd), ncol=size)
	W[sample.int(size^2, floor((1-density)*size^2))] <- 0
	W
}

test.rob.initenv <- function(W) {
	sapply(sd.test[["initenv"]], function(initenv.sd) 
		mean(robindex.initenv(W, a, dev.steps, dev.measure, initenv.sd, rep=rob.reps, log=log.robindex)))
}

test.rob.lateenv <- function(W) {
	sapply(sd.test[["lateenv"]], function(lateenv.sd) 
		mean(robindex.lateenv(W, a, dev.steps, dev.measure, lateenv.sd, rep=rob.reps, log=log.robindex)))
}

test.rob.initmut <- function(W) {
	sapply(sd.test[["initmut"]], function(initmut.sd) 
		mean(robindex.initmut(W, a, dev.steps, dev.measure, initmut.sd, rep=rob.reps, log=log.robindex)))
}

test.rob.latemut <- function(W) {
	sapply(sd.test[["latemut"]], function(latemut.sd) 
		mean(robindex.latemut(W, a, dev.steps, dev.measure, latemut.sd, rep=rob.reps, log=log.robindex)))
}

test.W <- function(W) {
	list(
		initenv = test.rob.initenv(W),
		lateenv = test.rob.lateenv(W),
		initmut = test.rob.initmut(W),
		latemut = test.rob.latemut(W)
	)
}

plotW <- function(testW, what="initenv", add=TRUE, xlim=NULL, ylim=NULL, type="l", col=col.phen[what], ...) {
	if (!add) {
		if (is.null(xlim)) xlim <- range(sd.test[[what]])
		if (is.null(ylim)) ylim <- range(testW[[what]])
		plot(NULL, xlim=xlim, ylim=ylim, xlab=as.expression(sigmas[what]), ylab=as.expression(phen[what]), log="x")
	}
	lines(sd.test[[what]], testW[[what]], col=col, type=type, ...)
}




pdf("figS6.pdf", width=6, height=12)
	layout(matrix(1:8, ncol=2))
	par(mar=c(4,4,3,0.5))

	# Random W matrices

	allW <- lapply(1:reps, function(i) randW(size=net.size, density=density, mean=reg.mean, sd=reg.sd))
	alltests <- mclapply(allW, test.W, mc.cores=mc.cores)
	
	for (tt in seq_along(alltests)) 
		plotW(alltests[[tt]], what="initenv", add=tt>1, col=tt, ylim=c(-40,-5))
	abline(v=default.initenv.sd, lty=3, col="darkgray")
	title("Random networks")
		
	for (tt in seq_along(alltests)) 
		plotW(alltests[[tt]], what="lateenv", add=tt>1, col=tt, ylim=c(-40,-3))		
	abline(v=default.lateenv.sd, lty=3, col="darkgray")
	
	for (tt in seq_along(alltests)) 
		plotW(alltests[[tt]], what="initmut", add=tt>1, col=tt, ylim=c(-40,-3))
	abline(v=default.initmut.sd, lty=3, col="darkgray")
		
	for (tt in seq_along(alltests)) 
		plotW(alltests[[tt]], what="latemut", add=tt>1, col=tt, ylim=c(-40,-3))
	abline(v=default.latemut.sd, lty=3, col="darkgray")


	# Evolved W matrices (from reference simulations, figG-null)
	
	pattern <- 'figG-null-\\d+.rds'
	gen <- NA
	
	allW <- lapply(
		sample(list.files(path="../cache", pattern=pattern, full.names=TRUE), reps, replace=FALSE), 
		function(ff) {
			ss <- readRDS(ff)
			if (is.na(gen) || !as.character(gen) %in% names(ss)) gen <- names(ss)[length(ss)]
			ss[[as.character(gen)]]$W
		})
	
	alltests <- mclapply(allW, test.W, mc.cores=mc.cores)

	for (tt in seq_along(alltests)) 
		plotW(alltests[[tt]], what="initenv", add=tt>1, col=tt, ylim=c(-40,-5))
	abline(v=default.initenv.sd, lty=3, col="darkgray")
	title("Evolved networks")
		
	for (tt in seq_along(alltests)) 
		plotW(alltests[[tt]], what="lateenv", add=tt>1, col=tt, ylim=c(-25,-3))		
	abline(v=default.lateenv.sd, lty=3, col="darkgray")
	
	for (tt in seq_along(alltests)) 
		plotW(alltests[[tt]], what="initmut", add=tt>1, col=tt, ylim=c(-25,-3))
	abline(v=default.initmut.sd, lty=3, col="darkgray")
		
	for (tt in seq_along(alltests)) 
		plotW(alltests[[tt]], what="latemut", add=tt>1, col=tt, ylim=c(-25,-3))
	abline(v=default.latemut.sd, lty=3, col="darkgray")

dev.off()
