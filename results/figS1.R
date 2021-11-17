#!/usr/bin/env Rscript

source("./terminology.R")
source("./defaults.R")

source("../src/randnetwork.R")
source("../src/robindex.R")

############### Options

param <- default

param$rob.reps      <- 1000

evolv.files.pattern <- 'figG-null-\\d+.rds'
evolv.gen <- NA

rep.show <- 8


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

##################### Functions


test.sensit.initenv <- function(W) {
	sapply(sd.test[["initenv"]], function(initenv.sd) 
		mean(robindex.initenv(W, param$a, param$dev.steps, param$dev.measure, initenv.sd, rep=param$rob.reps, log=param$log.robustness)))
}

test.sensit.lateenv <- function(W) {
	sapply(sd.test[["lateenv"]], function(lateenv.sd) 
		mean(robindex.lateenv(W, param$a, param$dev.steps, param$dev.measure, lateenv.sd, rep=param$rob.reps, log=param$log.robustness)))
}

test.sensit.initmut <- function(W) {
	sapply(sd.test[["initmut"]], function(initmut.sd) 
		mean(robindex.initmut(W, param$a, param$dev.steps, param$dev.measure, initmut.sd, rep=param$rob.reps, log=param$log.robustness)))
}

test.sensit.latemut <- function(W) {
	sapply(sd.test[["latemut"]], function(latemut.sd) 
		mean(robindex.latemut(W, param$a, param$dev.steps, param$dev.measure, latemut.sd, rep=param$rob.reps, log=param$log.robustness)))
}

test.sensit.W <- function(W) {
	list(
		initenv = test.sensit.initenv(W),
		lateenv = test.sensit.lateenv(W),
		initmut = test.sensit.initmut(W),
		latemut = test.sensit.latemut(W)
	)
}


plot.sensit.W <- function(testW, what="initenv", add=TRUE, xlim=NULL, ylim=NULL, ylab=as.expression(phen[what]), type="l", col=col.phen[what], ...) {
	if (!add) {
		if (is.null(xlim)) xlim <- range(sd.test[[what]])
		if (is.null(ylim)) ylim <- range(testW[[what]])
		plot(NULL, xlim=xlim, ylim=ylim, xlab=as.expression(sigmas[what]), ylab=ylab, log="x")
	}
	lines(sd.test[[what]], testW[[what]], col=col, type=type, ...)
}


######################## Figure

pdf("figS1.pdf", width=10/param$figscale, height=14/param$figscale, pointsize=param$pointsize)
	layout(matrix(1:8, ncol=2))
	par(mar=c(3,3,0.5,0.5), oma=c(0,2,2,0), cex=1, mgp=c(2,1,0))

	# Random W matrices

	allW.rand <- lapply(1:rep.show, function(i) randW(net.size=param$n, reg.mean=param$rand.mean, reg.sd=param$rand.sd, density=param$density))
	alltests.rand <- mclapply(allW.rand, test.sensit.W, mc.cores=param$mc.cores)
	
	for (tt in seq_along(alltests.rand)) 
		plot.sensit.W(alltests.rand[[tt]], what="initenv", add=tt>1, col=tt, ylim=c(-40,-5), ylab="")
	mtext(as.expression(phen["initenv"]), 2, col=default.cols["initenv"], line=2.5, xpd=NA)
	abline(v=param$initenv.sd, lty=3, col="darkgray")
	title("Random networks", xpd=NA, line=1)
		
	for (tt in seq_along(alltests.rand)) 
		plot.sensit.W(alltests.rand[[tt]], what="lateenv", add=tt>1, col=tt, ylim=c(-40,-3), ylab="")
	mtext(as.expression(phen["lateenv"]), 2, col=default.cols["lateenv"], line=2.5, xpd=NA)
	abline(v=param$lateenv.sd, lty=3, col="darkgray")
	
	for (tt in seq_along(alltests.rand)) 
		plot.sensit.W(alltests.rand[[tt]], what="initmut", add=tt>1, col=tt, ylim=c(-40,-3), ylab="")
	mtext(as.expression(phen["initmut"]), 2, col=default.cols["initmut"], line=2.5, xpd=NA)
	abline(v=param$initmut.sd, lty=3, col="darkgray")
		
	for (tt in seq_along(alltests.rand)) 
		plot.sensit.W(alltests.rand[[tt]], what="latemut", add=tt>1, col=tt, ylim=c(-40,-3), ylab="")
	mtext(as.expression(phen["latemut"]), 2, col=default.cols["latemut"], line=2.5, xpd=NA)
	abline(v=param$latemut.sd, lty=3, col="darkgray")


	# Evolved W matrices (from reference simulations, figG-null)
	

	allW.evolv <- lapply(
		sample(list.files(path=param$cache.dir, pattern=evolv.files.pattern, full.names=TRUE), rep.show, replace=FALSE), 
		function(ff) {
			ss <- readRDS(ff)
			if (is.na(evolv.gen) || !as.character(evolv.gen) %in% names(ss)) evolv.gen <- names(ss)[length(ss)]
			ss[[as.character(evolv.gen)]]$W
		})
	
	alltests.evolv <- mclapply(allW.evolv, test.sensit.W, mc.cores=param$mc.cores)

	for (tt in seq_along(alltests.evolv)) 
		plot.sensit.W(alltests.evolv[[tt]], what="initenv", add=tt>1, col=tt, ylim=c(-40,-5), ylab="")
	abline(v=param$initenv.sd, lty=3, col="darkgray")
	title("Evolved networks", xpd=NA, line=1)
		
	for (tt in seq_along(alltests.evolv)) 
		plot.sensit.W(alltests.evolv[[tt]], what="lateenv", add=tt>1, col=tt, ylim=c(-25,-3), ylab="")		
	abline(v=param$lateenv.sd, lty=3, col="darkgray")
	
	for (tt in seq_along(alltests.evolv)) 
		plot.sensit.W(alltests.evolv[[tt]], what="initmut", add=tt>1, col=tt, ylim=c(-25,-3), ylab="")
	abline(v=param$initmut.sd, lty=3, col="darkgray")
		
	for (tt in seq_along(alltests.evolv)) 
		plot.sensit.W(alltests.evolv[[tt]], what="latemut", add=tt>1, col=tt, ylim=c(-25,-3), ylab="")
	abline(v=param$latemut.sd, lty=3, col="darkgray")

dev.off()
