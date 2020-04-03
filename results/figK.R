#!/usr/bin/env Rscript

source("./commonpure.R")
source("./terminology.R")
source("./defaults.R")

library(parallel)
mc.cores <- default.mc.cores

use.cache <- TRUE

n.genes        <- default.n
sel.genes      <- default.nsel
s              <- c(rep(default.s, sel.genes), rep(0, n.genes-sel.genes))
W0             <- matrix(rnorm(n.genes^2, sd=default.initsd), ncol=n.genes)
reps           <- default.sim.reps
test.rep       <- default.rob.reps
N              <- default.N
G              <- 5000
every          <- round(G/100)
force.run      <- !use.cache

max.points     <- 20

phen <- c(
	fitness="Fitness",
    initenv=substitute(x~(y), list(x=TERM.ENVCAN.LONG, y=ABBRV.ENVCAN[[1]])),
    lateenv=substitute(x~(y), list(x=TERM.HOMEO.LONG, y=ABBRV.HOMEO[[1]])),
    initmut=substitute(x~(y), list(x=TERM.GENCAN.LONG, y=ABBRV.GENCAN[[1]])),
    latemut=substitute(x~(y), list(x=TERM.SOM.LONG, y=ABBRV.SOM[[1]])),
    stability=substitute(x~(y), list(x=TERM.STAB.LONG, y=ABBRV.STAB[[1]])))
col.phen <- c(fitness="black", initenv=COL.ENVCAN, lateenv=COL.HOMEO, initmut=COL.GENCAN, latemut=COL.SOM, stability=COL.STAB)
cols <- c('sim.ref'="black", 'sim.som'=COL.SOM, 'sim.iev'=COL.ENVCAN, 'sim.lev'=COL.HOMEO)
pchs <- c('sim.ref'=1, 'sim.som'=2, 'sim.iev'=0, 'sim.lev'=5)


allplots <- function(list.sim, what="fitness", xlim=NULL, ylim=NULL, xlab="Generations", ylab=what, lwd.highlight=NULL, FUN=mean, ...) {
	# Helper function to plot robustness time series
	
	if(is.null(xlim)) 
		xlim <- c(0, as.numeric(rev(names(list.sim[[1]]$mean))[1]))
	allv <- do.call(rbind, lapply(list.sim, function(x) do.call(rbind, lapply(x$mean, function(xx) FUN(xx[[what]])))))
	if(is.null(ylim)) 
		ylim <- c(min(allv), max(allv))
	
	
	plot(NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab="", main=as.expression(phen[what]), col.main=col.phen[what], ...)
	for (nss in names(list.sim)) {
		xpp <- seq_along(as.numeric(names(list.sim[[nss]]$mean)))
		xpp <- round(seq(xpp[1], xpp[length(xpp)], length.out=min(max.points, length(xpp))))
		lines(as.numeric(names(list.sim[[nss]]$mean))[xpp], sapply(list.sim[[nss]]$mean, function(x) FUN(x[[what]]))[xpp], col=cols[nss], lty=0, pch=pchs[nss], type="o", lwd=if(nss %in% lwd.highlight) 3 else 1)
	}
}

torun <- list(
	sim.ref = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, mut.rate=0.001),
		reps=reps, series.name="figK-ref", force.run=force.run),
	sim.som = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, mut.rate=0.001, som.mut.rate=0.1),
		reps=reps, series.name="figK-som", force.run=force.run),
	sim.iev = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, mut.rate=0.001, sim.initenv.sd=0.1),
		reps=reps, series.name="figK-iev", force.run=force.run),
	sim.lev = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, mut.rate=0.001, sim.lateenv.sd=0.1),
		reps=reps, series.name="figK-lev", force.run=force.run)
)

list.sim <- mclapply(torun, function(ff) ff(), mc.cores=min(length(torun), ceiling(mc.cores/reps)))


pdf("figK.pdf", width=10, height=12)
	layout(rbind(1:3,4:6))
	par(cex=1, mar=c(5, 2, 4, 1))
	allplots(list.sim, what="fitness")
	legend("topleft", pch=pchs, col=cols, legend=names(cols))
	allplots(list.sim, what="initenv", lwd.highlight="sim.iev")
	allplots(list.sim, what="lateenv", lwd.highlight="sim.lev")
	allplots(list.sim, what="initmut")
	allplots(list.sim, what="latemut", lwd.highlight="sim.som")
	allplots(list.sim, what="stability")
dev.off()
