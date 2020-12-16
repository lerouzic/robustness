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

pch.sim <- c(o=1, m=6, p=2, mm=10, pp=11, mp=9, pm=7)

plot.rel <- function(list.sim, what, G=NULL, pch=pch.sim[c("mm","mp","pm","pp")], xlab="Relative response", ...) {
	if (is.null(G)) G <- rev(names(list.sim[[1]]$mean))[1] # Default: last generation

	plot(NULL, xlim=c(-2,1), ylim=c(0.5,0.8 + length(what)), xlab=xlab, ylab="", yaxt="n", bty="n", ...)
	axis(2, at=seq_along(what), tick=FALSE, labels=as.expression(phen.expression[what]), las=2, mgp=c(3,0,0))
	abline(v=0, lty=3, col="darkgray")
	
	for (i in seq_along(what)) {
		ref.o <- mean(list.sim[["oo.o"]]$mean[[G]][[what[i]]])
		ref.m <- mean(list.sim[[paste0(default.shortcode[what[i]], ".m")]]$mean[[G]][[what[i]]])
		ref.p <- mean(list.sim[[paste0(default.shortcode[what[i]], ".p")]]$mean[[G]][[what[i]]])

		arrows(x0=(ref.m-ref.o)/(ref.p-ref.m), x1=(ref.p-ref.o)/(ref.p-ref.m), y0=i, code=3, angle=90, length=0.05, col=default.cols[what[i]], lwd=3)

		for (j in seq_along(what)) {
			if (j == i) next
			ss.name <- paste0(default.shortcode[what[i]], ".", default.shortcode[what[j]])
			ss <- sapply(names(pch), function(nn) mean(list.sim[[paste0(ss.name, ".", nn)]]$mean[[G]][[what[i]]]))
			points((ss-ref.o)/(ref.p-ref.m), rep(i+0.1*j, length(ss)), pch=pch[names(ss)], col=default.cols[what[j]])
		}
	}
	legend("topright", pch=pch, legend=c("target - & corr -     ", "target - & corr +     ", "target + & corr -", "target + & corr +    "), horiz=TRUE, bty="n", xpd=NA, cex=0.8)
}

plot.ts <- function(list.sim, what, w1, w2, xlab="Generation", ylab=as.expression(phen.expression[what]), ylim=NULL, lwd=2, ...) {
	Gmax <- as.numeric(rev(names(list.sim[[1]]$mean))[1])
	xx <- unique(round(seq(from=1, to=length(list.sim[[1]]$mean), length.out=15)))
	
	sw1 <- default.shortcode[w1]
	sw2 <- default.shortcode[w2]
	
	ref.o <- sapply(list.sim[["oo.o"]]$mean, function(x) mean(x[[what]]))
	ref.m <- sapply(list.sim[[paste0(sw1, ".m")]]$mean, function(x) mean(x[[what]]))
	ref.p <- sapply(list.sim[[paste0(sw1, ".p")]]$mean, function(x) mean(x[[what]]))
	
	sim.mm <- sapply(list.sim[[paste0(sw1, ".", sw2, ".mm")]]$mean, function(x) mean(x[[what]]))
	sim.mp <- sapply(list.sim[[paste0(sw1, ".", sw2, ".mp")]]$mean, function(x) mean(x[[what]]))
	sim.pm <- sapply(list.sim[[paste0(sw1, ".", sw2, ".pm")]]$mean, function(x) mean(x[[what]]))
	sim.pp <- sapply(list.sim[[paste0(sw1, ".", sw2, ".pp")]]$mean, function(x) mean(x[[what]]))
	
	if (is.null(ylim))
		ylim <- range(c(ref.m, ref.p, sim.mm, sim.pp))
	
	plot(NULL, xlim=c(1, Gmax), ylim=ylim, xlab=xlab, ylab=ylab, ...)
	lines(as.numeric(names(ref.o))[xx], ref.o[xx], pch=pch.sim["o"], lwd=lwd)
	lines(as.numeric(names(ref.p))[xx], ref.p[xx], lty=1, col=default.cols[w1], lwd=lwd) #pch=pch.sim["p"])
	lines(as.numeric(names(ref.m))[xx], ref.m[xx], lty=1, col=default.cols[w1], lwd=lwd) #pch=pch.sim["m"])
	
	points(as.numeric(names(sim.mm))[xx], sim.mm[xx], pch=pch.sim["mm"], col=default.cols[w2])
	points(as.numeric(names(sim.mp))[xx], sim.mp[xx], pch=pch.sim["mp"], col=default.cols[w2])
	points(as.numeric(names(sim.pm))[xx], sim.pm[xx], pch=pch.sim["pm"], col=default.cols[w2])
	points(as.numeric(names(sim.pp))[xx], sim.pp[xx], pch=pch.sim["pp"], col=default.cols[w2])

}

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
	
	rat <- sapply(ans, "[", "ratio")
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

# It is more convenient to have all combinations in the list.sim variable (should not take more memory)
for (i in 1:(length(indwhat)-1))
	for (j in (i+1):length(indwhat)) {
		nc <- paste0(names(indwhat)[i], ".", names(indwhat)[j], ".")
		inc <- paste0(names(indwhat)[j], ".", names(indwhat)[i], ".")
		list.sim[[paste0(inc, "pp")]] <- list.sim[[paste0(nc, "pp")]]
		list.sim[[paste0(inc, "mp")]] <- list.sim[[paste0(nc, "pm")]]
		list.sim[[paste0(inc, "pm")]] <- list.sim[[paste0(nc, "mp")]]
		list.sim[[paste0(inc, "mm")]] <- list.sim[[paste0(nc, "mm")]]
	}

# Realized evolvabilities
real.evolv <- real.evolv.sd <- matrix(NA, ncol=length(indwhat)-1, nrow=length(indwhat)-1)
colnames(real.evolv) <- colnames(real.evolv.sd) <- names(indwhat)[1:(length(indwhat)-1)]
rownames(real.evolv) <- rownames(real.evolv.sd) <- names(indwhat)[2:length(indwhat)]

for (i1 in 1:(length(indwhat)-1))
	for (i2 in(i1+1):length(indwhat)) {
		n1 <- names(indwhat)[i1]
		n2 <- names(indwhat)[i2]
		rr <- rel.response(list.sim, n1, n2, indwhat[i1], indwhat[i2], normalize=FALSE)
		real.evolv[n2, n1] <- rr$mean.ratio
		real.evolv.sd[n2, n1] <- rr$sd.ratio
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


pdf("fig4.pdf", width=8, height=5)
	layout(cbind(c(1,1,1),c(2,3,4)), widths=c(2,1))
	par(mar=c(4, 10, 1, 1), cex=0.75, oma=c(2,0,0,0))
	plot.rel(list.sim, indwhat, xlab="")
	mtext(side=1, "Relative response", xpd=NA, line=3, cex=0.75)
	
	arrows(x0=0.5, y0=5.4, x1=1.2, y1=5.4, lwd=2, xpd=NA, length=0.1, col=default.cols["latemut"])
	arrows(x0=0.5, y0=3.2, x1=1.3, y1=3.25, lwd=2, xpd=NA, length=0.1, col=default.cols["lateenv"])
	arrows(x0=0.7, y0=1.5, x1=1.3, y1=1.3, lwd=2, xpd=NA, length=0.1, col=default.cols["stability"])
	
	par(mar=c(2, 4, 1, 1))
	plot.ts(list.sim, "stability", "stability", "latemut", xlab="", mgp=c(2,1,0))
	plot.ts(list.sim, "initmut", "initmut", "lateenv", xlab="", mgp=c(2,1,0), ylim=c(-11,-6))
	plot.ts(list.sim, "initenv", "initenv", "stability", xlab="", mgp=c(2,1,0))
	mtext(side=1, "Generations", xpd=NA, line=2.5, cex=0.75)
	
dev.off()



