#!/usr/bin/env Rscript

source("./commonsim.R")
source("./terminology.R")
source("./defaults.R")
source("../src/robindex.R")
cache.dir <- "../cache"

library(parallel)
library(abind)
library(ellipse)

mc.cores <- default.mc.cores

use.cache <- TRUE

n.genes           <- default.n
sel.genes         <- default.nsel    
s                 <- c(rep(default.s, sel.genes), rep(0, n.genes-sel.genes))
W0                <- NA
reps              <- default.sim.reps
test.rep          <- default.rob.reps
grad.effect       <- 0.01
N                 <- default.N 
n                 <- default.n
a                 <- default.a
G                 <- 10000
every             <- round(G/100)
force.run         <- !use.cache

M.reps            <- 100
M.nbmut           <- 5
gen.display       <- 1000

col.sim  <- c(oo="black", ie=COL.ENVCAN, le=COL.HOMEO, im=COL.GENCAN, lm=COL.SOM, st=COL.STAB)
pch.sim <- c(o=1, m=6, p=2, mm=10, pp=11, mp=9, pm=7)
xylim <- list(
	initenv = c(-38, -12),
	lateenv = c(- 8, - 5.8),
	initmut = c(-8.5,- 6),
	latemut = c(-10,- 7), 
	stability=c(-38, -12))

list.mean <- function(ll) {
	arr <- do.call(abind, c(ll, list(along=3)))
	rowMeans(arr, dims=2)
}


################# 
plot2traits <- function(list.sim, list.M=NULL, what.x, what.y, xlim=NULL, ylim=NULL, xlab=as.expression(phen.expression[what.x]), ylab=as.expression(phen.expression[what.y]), FUN=mean, mut.rate=1, ...) {
	if (is.null(xlim))
		xlim <- range(do.call(rbind, lapply(list.sim, function(x) do.call(rbind, lapply(x$mean, function(xx) FUN(xx[[what.x]]))))))
	if (is.null(ylim))
		ylim <- range(do.call(rbind, lapply(list.sim, function(x) do.call(rbind, lapply(x$mean, function(xx) FUN(xx[[what.y]]))))))
	
	plot(NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
	
	ref.x <- sapply(list.sim[["oo.o"]]$mean, function(x) FUN(x[[what.x]]))[1]
	ref.y <- sapply(list.sim[["oo.o"]]$mean, function(x) FUN(x[[what.y]]))[1]
	if (!is.null(list.M)) {
		lines(ellipse(mut.rate*list.M[["oo.o"]]$mean.M[c(what.x, what.y), c(what.x, what.y)], centre=c(ref.x, ref.y)), lty=2)
		lines(ellipse(mut.rate*list.M[["oo.o"]]$mean.Mcond[c(what.x, what.y), c(what.x, what.y)], centre=c(ref.x, ref.y)))
	}
		
	for (nss in names(list.sim)) {
		nnss <- strsplit(nss, split="\\.")[[1]]
		pch <- col <- NA
		if (length(nnss) == 2) {
			if (nnss[1] == default.shortcode[what.x])
				col <- col.sim[default.shortcode[what.x]]
			if (nnss[1] == default.shortcode[what.y])
				col <- col.sim[default.shortcode[what.y]]
			pch <- pch.sim[nnss[2]]
		}
		if (length(nnss) == 3) {
			if (nnss[1] == default.shortcode[what.x] && nnss[2] == default.shortcode[what.y]) {
				col <- "darkgray"
				pch <- pch.sim[nnss[3]]
			}
		}

		if (is.na(pch) || is.na(col)) next
		gen.numbers <- as.numeric(names(list.sim[[nss]]$mean))
		xpp <- seq_along(gen.numbers[gen.numbers <= gen.display])
#~ 		xpp <- round(seq(xpp[1], xpp[length(xpp)], length.out=min(max.points, length(xpp))))
		lines(
			x = c(ref.x, sapply(list.sim[[nss]]$mean, function(x) FUN(x[[what.x]]))[xpp]),
			y = c(ref.y, sapply(list.sim[[nss]]$mean, function(x) FUN(x[[what.y]]))[xpp]),
			col = col)
		points(x=ref.x, y=ref.y, pch=pch.sim["o"])
		points(
			x = sapply(list.sim[[nss]]$mean, function(x) FUN(x[[what.x]]))[xpp[length(xpp)]],
			y = sapply(list.sim[[nss]]$mean, function(x) FUN(x[[what.y]]))[xpp[length(xpp)]],
			col = col, 
			pch = pch,
			cex = 2)
	}
}

plot2traits.legend <- function(col.x, col.y, corr.M=0.5, evol.dist=1.5, ...) {
	plot(NULL, xlim=c(-2, 2), ylim=c(-2, 2), xaxt="n", yaxt="n", xlab="", ylab="", bty="n", ...)
	
	lines(ellipse(0.2*rbind(c(1, corr.M), c(corr.M, 1))), lty=2)
	lines(ellipse(0.15*rbind(c(1, corr.M), c(corr.M, 1))),)
	legend("bottom", lty=c(2,1), legend=c("M (x10)", "Mcond (x10)"), horiz=TRUE, bty="n")
	points(0, 0, pch=pch.sim["o"])
	points(x=c(0,0), y=evol.dist*c(-1,1), pch=pch.sim[c("m","p")], col=col.y, cex=2)
	arrows(x0=c(0,0), y0=c(0,0), x1=c(0,0), y1=0.8*evol.dist*c(-1,1), col=col.y, length=0.1)
	points(x=evol.dist*c(-1,1), y=c(0,0), pch=pch.sim[c("m","p")], col=col.x, cex=2)
	arrows(x0=c(0,0), y0=c(0,0), x1=0.8*evol.dist*c(-1,1), y1=c(0,0), col=col.x, length=0.1)
	points(x=evol.dist*c(-1,-1,1,1), y=evol.dist*c(-1,1,1,-1), pch=pch.sim[c("mm","mp","pp","pm")], col="darkgray", cex=2)
	arrows(x0=c(0,0,0,0), y0=c(0,0,0,0), x1=0.8*evol.dist*c(-1,-1,1,1), y1=0.8*evol.dist*c(-1,1,1,-1), col="darkgray",length=0.1)
}


# Computes gradient verctors for one and two selected traits
gradvec1 <- function(i, grd) c(rep(0, i-1), grd, rep(0, length(default.shortcode)-i))
gradvec2 <- function(i1, i2, grd1, grd2) c(rep(0, i1-1), grd1, rep(0, i2-i1-1), grd2, rep(0, length(default.shortcode)-i2))

# Make a list of simulation names and gradients
torun <- list(
	oo.o = list(series.name="figG-null", grad.rob=c(0,0,0,0,0)))
for (i in seq_along(default.shortcode)) {
	nm <- default.shortcode[i]
	torun[[paste(nm, "m", sep=".")]] <- list(series.name=paste("figG", nm, "m", sep="-"), grad.rob=gradvec1(i, -grad.effect))
	torun[[paste(nm, "p", sep=".")]] <- list(series.name=paste("figG", nm, "p", sep="-"), grad.rob=gradvec1(i,  grad.effect))
}
for (i1 in 1:(length(default.shortcode)-1))
	for (i2 in (i1+1):length(default.shortcode)) {
		nm1 <- default.shortcode[i1]
		nm2 <- default.shortcode[i2]
		torun[[paste(nm1, nm2, "mm", sep=".")]] <- list(series.name=paste("figI", nm1, nm2, "mm", sep="-"), grad.rob=gradvec2(i1, i2, -grad.effect, -grad.effect))
		torun[[paste(nm1, nm2, "mp", sep=".")]] <- list(series.name=paste("figI", nm1, nm2, "mp", sep="-"), grad.rob=gradvec2(i1, i2, -grad.effect,  grad.effect))
		torun[[paste(nm1, nm2, "pm", sep=".")]] <- list(series.name=paste("figI", nm1, nm2, "pm", sep="-"), grad.rob=gradvec2(i1, i2,  grad.effect, -grad.effect))
		torun[[paste(nm1, nm2, "pp", sep=".")]] <- list(series.name=paste("figI", nm1, nm2, "pp", sep="-"), grad.rob=gradvec2(i1, i2,  grad.effect,  grad.effect))
	}

list.sim <- mclapply(torun, function(ff) 
	sim.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=ff$grad.rob), reps=reps, series.name=ff$series.name, force.run=force.run, mc.cores=reps), 
	mc.cores=ceiling(mc.cores/reps))
	
# Computing M matrices (so far, only need the control for the first generation)
list.M <- mclapply(list.sim["oo.o"], function(sim.series) {
		W0.all <- lapply(sim.series$full, function(x) x[[1]]$W)
		M0.all <- mclapply(W0.all, robindex.Mmatrix, a=a, dev.steps=default.dev.steps, mut.sd=default.sim.mutsd, mut.correlated=default.mut.correlated, test.initmut.sd=default.initmut.sd, test.latemut.sd=default.latemut.sd, nbmut=M.nbmut, test.initenv.sd=default.initenv.sd, test.lateenv.sd=default.lateenv.sd, test.rep=test.rep, rep=M.reps, log.robustness=default.log.robustness, include.expr=TRUE, mc.cores=mc.cores)
		list(
			full      = M0.all, 
			mean.M    = list.mean(lapply(M0.all, function(m) m$vcov[names(default.shortcode), names(default.shortcode)])),
			mean.Mcond= list.mean(lapply(M0.all, function(m) conditional(m$vcov, paste0("expr", seq_len(sel.genes)))[names(default.shortcode), names(default.shortcode)]))) 
	}, mc.cores=1)
		

# It is more convenient to have all combinations in the list.sim variable (should not take more memory)
for (i in 1:(length(default.shortcode)-1))
	for (j in (i+1):length(default.shortcode)) {
		nc <- paste0(default.shortcode[i], ".", default.shortcode[j], ".")
		inc <- paste0(default.shortcode[j], ".", default.shortcode[i], ".")
		list.sim[[paste0(inc, "pp")]] <- list.sim[[paste0(nc, "pp")]]
		list.sim[[paste0(inc, "mp")]] <- list.sim[[paste0(nc, "pm")]]
		list.sim[[paste0(inc, "pm")]] <- list.sim[[paste0(nc, "mp")]]
		list.sim[[paste0(inc, "mm")]] <- list.sim[[paste0(nc, "mm")]]
	}

pdf("fig5.pdf", width=10, height=10)
	lm <- matrix(0, ncol=4, nrow=4)
	lm[lower.tri(lm, diag=TRUE)] <- 1:10
	lm[1,4] <- 11
	layout(lm)
	
	par(mar=c(0.5, 0.5, 0.5, 0.5), oma=c(4, 4, 0, 2))
	
	for (i in 1:4) {
		for (j in (i+1):5) {
			plot2traits(list.sim, list.M, names(default.shortcode)[i], names(default.shortcode)[j], xaxt="n", yaxt="n", xlab="", ylab="", xlim=xylim[[i]], ylim=xylim[[j]])
			if (i==1) mtext(as.expression(phen.expression[j]), 2, line=2, xpd=NA, col=default.cols[j])
			if (j==i+1) axis(4)
			if (j==5) { axis(1); mtext(as.expression(phen.expression[i]), 1, line=3, xpd=NA, col=default.cols[i]) }
		}
	}
	plot2traits.legend(default.cols[4], default.cols[2])
	mtext("Direction of selection", 1, cex=1.5, line=1)
dev.off()
