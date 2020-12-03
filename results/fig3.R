#!/usr/bin/env Rscript

source("./commonpure.R")
source("./terminology.R")
source("./defaults.R")

library(parallel)
mc.cores <- default.mc.cores

use.cache <- TRUE

n.genes          <- default.n
sel.genes        <- default.nsel
s                <- c(rep(default.s, sel.genes), rep(0, n.genes-sel.genes))
W0               <- NA # Random W0

reps             <- default.sim.reps
test.rep         <- default.rob.reps
grad.effect      <- 0.01
N                <- default.N
G                <- 10000
every            <- max(1, round(G/100))
force.run        <- !use.cache

phen <- c(list(fitness="Fitness"), phen.expression)
col.sim <- c(oo="black", ie=COL.ENVCAN, le=COL.HOMEO, im=COL.GENCAN, lm=COL.SOM, st=COL.STAB)
col.phen <- c(fitness="black", initenv=COL.ENVCAN, lateenv=COL.HOMEO, initmut=COL.GENCAN, latemut=COL.SOM, stability=COL.STAB)
lty.sim <- c(p=0, m=0, o=0)
pch.sim <- c(p=2, m=6, o=1)
max.points <- 20

allplots <- function(list.sim, what="fitness", xlim=NULL, ylim=NULL, xlab="Generations", ylab=what, lwd=3, FUN=mean, focal="oo", ...) {
	# Helper function to plot robustness time series
	
	if(is.null(xlim)) 
		xlim <- c(0, as.numeric(rev(names(list.sim[[1]]$mean))[1]))
	allv <- do.call(rbind, lapply(list.sim, function(x) do.call(rbind, lapply(x$mean, function(xx) FUN(xx[[what]])))))
	if(is.null(ylim)) 
		ylim <- c(min(allv), max(allv))
	
	# Reordering the data so that the "focal" traits are plotted at the end and will overlap the others
	if (length(focal) > 0) list.sim <- list.sim[order(sapply(strsplit(names(list.sim), split="\\."), function(x) x[1]) %in% focal)]
	
	plot(NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab="", main=as.expression(phen[what]), col.main=col.phen[what], ...)
	for (nss in names(list.sim)) {
		nnss <- strsplit(nss, split="\\.")[[1]]
		xpp <- seq_along(as.numeric(names(list.sim[[nss]]$mean)))
		xpp <- round(seq(xpp[1], xpp[length(xpp)], length.out=min(max.points, length(xpp))))
		if(any(is.na(sapply(list.sim[[nss]]$mean, function(x) FUN(x[[what]]))[xpp]))) browser()
		lines(as.numeric(names(list.sim[[nss]]$mean))[xpp], sapply(list.sim[[nss]]$mean, function(x) FUN(x[[what]]))[xpp], col=col.sim[nnss[1]], lty=lty.sim[nnss[2]], pch=pch.sim[nnss[2]], type="o", lwd=if(nnss[1] %in% focal) 3 else 1)
	}
}

torun <- list(
	oo.o = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,0,0,0,0)), 
		reps=reps, series.name="figG-null", force.run=force.run), 
	ie.m = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(-grad.effect,0,0,0,0)), 
		reps=reps, series.name="figG-ie-m", force.run=force.run),
	le.m = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,-grad.effect,0,0,0)), 
		reps=reps, series.name="figG-le-m", force.run=force.run),
	im.m = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,0,-grad.effect,0,0)), 
		reps=reps, series.name="figG-im-m", force.run=force.run),
	lm.m = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,0,0,-grad.effect,0)), 
		reps=reps, series.name="figG-lm-m", force.run=force.run),
	st.m = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,0,0,0,-grad.effect)), 
		reps=reps, series.name="figG-st-m", force.run=force.run),
	ie.p = function()pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(grad.effect,0,0,0,0)), 
		reps=reps, series.name="figG-ie-p", force.run=force.run),
	le.p = function()pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,grad.effect,0,0,0)), 
		reps=reps, series.name="figG-le-p", force.run=force.run),
	im.p = function()pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,0,grad.effect,0,0)), 
		reps=reps, series.name="figG-im-p", force.run=force.run),
	lm.p = function()pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,0,0,grad.effect,0)), 
		reps=reps, series.name="figG-lm-p", force.run=force.run),
	st.p = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,0,0,0,grad.effect)), 
		reps=reps, series.name="figG-st-p", force.run=force.run)
)

list.sim <- mclapply(torun, function(ff) ff(), mc.cores=min(length(torun), ceiling(mc.cores/reps)))


pdf("fig3.pdf", width=12, height=8)
	layout(rbind(1:3,4:6))
	par(cex=1, mar=c(5, 2, 4, 1))
	allplots(list.sim, what="initenv", focal=c("oo","ie"))
	allplots(list.sim, what="lateenv", focal=c("oo","le"))
	allplots(list.sim, what="initmut", focal=c("oo","im"))
	allplots(list.sim, what="latemut", focal=c("oo","lm"))
	allplots(list.sim, what="stability", focal=c("oo","st"))
dev.off()
