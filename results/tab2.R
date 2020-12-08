#!/usr/bin/env Rscript

source("./commonsim.R")
source("./terminology.R")
source("./defaults.R")
source("../src/robindex.R")
cache.dir <- "../cache"

library(parallel)
library(abind)

mc.cores <- default.mc.cores

use.cache <- TRUE

n.genes           <- default.n
sel.genes         <- default.nsel    
s                 <- c(rep(default.s, sel.genes), rep(0, n.genes-sel.genes))
W0                <- NA
reps              <- 10
test.rep          <- default.rob.reps
mut.rep           <- 1000
grad.effect       <- 0.01
rand.mean         <- default.rand.mean
rand.sd           <- default.rand.sd.sim
N                 <- default.N 
n                 <- default.n
a                 <- default.a
sim.mutsd         <- default.sim.mutsd
dev.steps         <- default.dev.steps
dev.measure       <- default.dev.measure
initenv.sd        <- default.initenv.sd
lateenv.sd        <- default.lateenv.sd
initmut.sd        <- default.initmut.sd
latemut.sd        <- default.latemut.sd
log.robustness    <- default.log.robustness
mut.correlated    <- default.mut.correlated
G                 <- 10000
every             <- round(G/100)
force.run         <- !use.cache

rel.response.points <- function(allpoints, normalize=FALSE) {
	stopifnot(nrow(allpoints)==4, ncol(allpoints)==2)
	
	if (normalize) {
		allpoints[,1] <- allpoints[,1]-allpoints["mm",1]
		allpoints[,2] <- allpoints[,2]-allpoints["mm",2]
		allpoints[,1] <- allpoints[,1] / allpoints["pp",1]
		allpoints[,2] <- allpoints[,2] / allpoints["pp",2]
	}
	
	# Direction of the most evolvability: mm to pp
	a <- (allpoints["pp",2] - allpoints["mm",2])/(allpoints["pp",1] - allpoints["mm",1])
	ap <- - 1/a # perpendicular direction
	
	mpp.x <- (allpoints["mp",1] + ap*allpoints["mp",2])/(1+ap^2)
	mpp.y <- ap*(allpoints["mp",1] + ap*allpoints["mp",2])/(1+ap^2)
	
	pmp.x <- (allpoints["pm",1] + ap*allpoints["pm",2])/(1+ap^2)
	pmp.y <- ap*(allpoints["pm",1] + ap*allpoints["pm",2])/(1+ap^2)
	
	dist.mm.pp <- sqrt((allpoints["mm",1]-allpoints["pp",1])^2 + (allpoints["mm",2] - allpoints["pp",2])^2)
	dist.mp.pm <- sqrt((mpp.x - pmp.x)^2 + (mpp.y - pmp.y)^2)
	return(c(max.dist=unname(dist.mm.pp), min.dst=unname(dist.mp.pm), ratio=unname(dist.mp.pm/dist.mm.pp)))
}

rel.response.onesim <- function(mm, mp, pm, pp, what1, what2, G=NULL, normalize=FALSE) {
	if (is.null(G)) G <- rev(names(mm))[1] # Default: last generation
	
	allpoints <- rbind( mm = c(mean(mm[[G]][[what1]]), mean(mm[[G]][[what2]])), 
						mp = c(mean(mp[[G]][[what1]]), mean(mp[[G]][[what2]])), 
						pm = c(mean(pm[[G]][[what1]]), mean(pm[[G]][[what2]])), 
						pp = c(mean(pp[[G]][[what1]]), mean(pp[[G]][[what2]])))
						
	rel.response.points(allpoints)
}

rel.response <- function(list.sim, r1, r2, what1, what2, G=NULL, normalize=FALSE) {
	
	ssim <- c(paste(r1, r2, c("mm","mp","pm","pp"), sep="."))
		
	if (any(!ssim %in% names(list.sim))) return(list(mean.ratio=NA, sd.ratio=NA, se.ratio=NA))
	
	# loop over replicates
	rep <- length(list.sim[[ssim[1]]]$full)
	ans <- lapply(1:rep, function(i) 
		rel.response.onesim(mm = list.sim[[ssim[1]]]$full[[i]], 
		                    mp = list.sim[[ssim[2]]]$full[[i]], 
		                    pm = list.sim[[ssim[3]]]$full[[i]], 
		                    pp = list.sim[[ssim[4]]]$full[[i]],
		                    what1=what1, what2=what2, G=G, normalize=normalize))
	
	rat <- sapply(ans, "[", "ratio.mp")
	list(mean.ratio = mean(rat), sd.ratio = sd(rat), se.ratio = sd(rat)/sqrt(rep))
}

mutate <- function(W, mut.sd, nb.mut=1) {
	for (nn in 1:nb.mut) {
		which.mut <- sample(size=1,  which(W != 0)) # Bug if only one W != 0
		W[which.mut] <- rnorm(1, mean=if(mut.correlated) W[which.mut] else 0, sd=mut.sd)
	}
	W
}

randW0 <- function(theta=runif(n)) { # Make a W matrix with the same algorithm as in the simulations
	# NA for W0 is replaced by a random matrix at the correct equilibrium
	WW <- matrix(rnorm(length(theta)^2, rand.mean, rand.sd), ncol=length(theta))
	for (i in 1:nrow(WW)) WW[i,sample(1:ncol(WW), 1)] <- NA
	targetW(WW, target=theta, a=a)
}

M.mat <- function(W, nbmut=2, mc.cores=1, mut.sd=sim.mutsd) {
	tt <- mclapply(1:mut.rep, function(rr) {
		myW <- mutate(W, mut.sd=mut.sd)
		c(	initenv  =mean(robindex.initenv  (myW, a, dev.steps, dev.measure, initenv.sd, rep=test.rep, log=log.robustness)),
			lateenv  =mean(robindex.lateenv  (myW, a, dev.steps, dev.measure, lateenv.sd, rep=test.rep, log=log.robustness)),
			initmut  =mean(robindex.initmut  (myW, a, dev.steps, dev.measure, initmut.sd, rep=test.rep, log=log.robustness)),
			latemut  =mean(robindex.latemut  (myW, a, dev.steps, dev.measure, latemut.sd, rep=test.rep, log=log.robustness)),
			stability=mean(robindex.stability(myW, a, dev.steps, dev.measure, log=log.robustness)))
		}, mc.cores=mc.cores)
	var(do.call(rbind, tt))
}

avg.M.mat <- function(W.rep = 100) {
	ans <- mclapply(1:W.rep, function(rr) {
		W <- randW0()
		M.mat(W)
	}, mc.cores=mc.cores)
	arr <- do.call(abind, c(ans, list(along=3)))
	rowMeans(arr, dims=2)
}


indx <- c(ie=1, le=2, im=3, lm=4, st=5)
gradvec1 <- function(i, grd) c(rep(0, i-1), grd, rep(0, length(indx)-i))
gradvec2 <- function(i1, i2, grd1, grd2) c(rep(0, i1-1), grd1, rep(0, i2-i1-1), grd2, rep(0, length(indx)-i2))

# Make a list of simulation names and gradients

torun <- list(
	oo.o = list(series.name="figG-null", grad.rob=c(0,0,0,0,0)))
for (i in indx) {
	nm <- names(indx)[i]
	torun[[paste(nm, "m", sep=".")]] <- list(series.name=paste("figG", nm, "m", sep="-"), grad.rob=gradvec1(i, -grad.effect))
	torun[[paste(nm, "p", sep=".")]] <- list(series.name=paste("figG", nm, "p", sep="-"), grad.rob=gradvec1(i,  grad.effect))
}
for (i1 in 1:(length(indx)-1))
	for (i2 in (i1+1):length(indx)) {
		nm1 <- names(indx)[i1]
		nm2 <- names(indx)[i2]
		torun[[paste(nm1, nm2, "mm", sep=".")]] <- list(series.name=paste("figI", nm1, nm2, "mm", sep="-"), grad.rob=gradvec2(i1, i2, -grad.effect, -grad.effect))
		torun[[paste(nm1, nm2, "mp", sep=".")]] <- list(series.name=paste("figI", nm1, nm2, "mp", sep="-"), grad.rob=gradvec2(i1, i2, -grad.effect,  grad.effect))
		torun[[paste(nm1, nm2, "pm", sep=".")]] <- list(series.name=paste("figI", nm1, nm2, "pm", sep="-"), grad.rob=gradvec2(i1, i2,  grad.effect, -grad.effect))
		torun[[paste(nm1, nm2, "pp", sep=".")]] <- list(series.name=paste("figI", nm1, nm2, "pp", sep="-"), grad.rob=gradvec2(i1, i2,  grad.effect,  grad.effect))
	}

list.sim <- mclapply(torun, function(ff) 
	sim.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=ff$grad.rob), reps=reps, series.name=ff$series.name, force.run=force.run, mc.cores=reps), 
	mc.cores=ceiling(mc.cores/reps))


indwhat <- c(ie="initenv", le="lateenv", im="initmut", lm="latemut", st="stability")

# Realized evolvabilities
real.evolv <- matrix("", ncol=length(indwhat)-1, nrow=length(indwhat)-1)
colnames(real.evolv) <- names(indwhat)[1:(length(indwhat)-1)]
rownames(real.evolv) <- names(indwhat)[2:length(indwhat)]

for (i1 in 1:(length(indwhat)-1))
	for (i2 in(i1+1):length(indwhat)) {
		n1 <- names(indwhat)[i1]
		n2 <- names(indwhat)[i2]
		rr <- rel.response(list.sim, n1, n2, indwhat[i1], indwhat[i2], normalize=FALSE)
		real.evolv[n2, n1] <- paste0(round(rr$mean.ratio, digits=2), " +/- ", round(rr$sd.ratio, digits=2))
	}


# Mutational evolvabilities
mut.evolv <- matrix(NA, ncol=length(indwhat)-1, nrow=length(indwhat)-1)
colnames(mut.evolv) <- names(indwhat)[1:(length(indwhat)-1)]
rownames(mut.evolv) <- names(indwhat)[2:length(indwhat)]

mut.cache.file <- file.path(cache.dir, "figI-mut.rds")
if (use.cache && file.exists(mut.cache.file)) {
	avgM <- readRDS(mut.cache.file)
} else {
	avgM <- avg.M.mat(W.rep=100) # average from 100 starting W0 matrices
	saveRDS(avgM, mut.cache.file)
}

for (i1 in 1:(length(indwhat)-1))
	for (i2 in(i1+1):length(indwhat)) {
		resp <- rbind(
			mm = (avgM %*% gradvec2(i1, i2, -1, -1))[c(i1,i2)],
			mp = (avgM %*% gradvec2(i1, i2, -1,  1))[c(i1,i2)],
			pm = (avgM %*% gradvec2(i1, i2,  1, -1))[c(i1,i2)],
			pp = (avgM %*% gradvec2(i1, i2,  1,  1))[c(i1,i2)])
		n1 <- names(indwhat)[i1]
		n2 <- names(indwhat)[i2]
		rr <- rel.response.points(resp, normalize=FALSE)
		mut.evolv[n2, n1] <- rr["ratio"] 
	}

format(mut.evolv, nsmall=2, digits=2)
