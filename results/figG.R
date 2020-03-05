#!/usr/bin/env Rscript

source("./commonpure.R")
source("./terminology.R")

n.genes <- 6
sel.genes <- 3
s <- c(rep(10, sel.genes), rep(0, n.genes-sel.genes))
W0 <- matrix(rnorm(n.genes^2, sd=0.000001), ncol=n.genes)
reps <- 3
test.rep <- 10
grad.effect <- 0.01
N <- 500
G <- 1000
force.run <- FALSE

col.sim <- c(oo="black", ie=COL.ENVCAN, le=COL.HOMEO, im=COL.GENCAN, lm=COL.SOM, st=COL.STAB)
col.phen <- c(fitness="black", initenv=COL.ENVCAN, lateenv=COL.HOMEO, initmut=COL.GENCAN, latemut=COL.SOM, stability=COL.STAB)
lty.sim <- c(p=0, m=0, o=0)
pch.sim <- c(p=2, m=6, o=1)
max.points <- 20
phen <- c(
	fitness="Fitness",
    initenv=substitute(x~(y), list(x=TERM.ENVCAN.LONG, y=ABBRV.ENVCAN[[1]])),
    lateenv=substitute(x~(y), list(x=TERM.HOMEO.LONG, y=ABBRV.HOMEO[[1]])),
    initmut=substitute(x~(y), list(x=TERM.GENCAN.LONG, y=ABBRV.GENCAN[[1]])),
    latemut=substitute(x~(y), list(x=TERM.SOM.LONG, y=ABBRV.SOM[[1]])),
    stability=substitute(x~(y), list(x=TERM.STAB.LONG, y=ABBRV.STAB[[1]])))

allplots <- function(list.sim, what="fitness", xlim=NULL, ylim=NULL, xlab="Generations", ylab=what, lwd=3, FUN=mean, ...) {
	if(is.null(xlim)) xlim <- c(0, as.numeric(rev(names(list.sim[[1]]$mean))[1]))
	allv <- do.call(rbind, lapply(list.sim, function(x) do.call(rbind, lapply(x$mean, function(xx) FUN(xx[[what]])))))
	if(is.null(ylim)) ylim <- c(min(allv), max(allv))
	plot(NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab="", main=as.expression(phen[what]), col.main=col.phen[what], ...)
	for (nss in names(list.sim)) {
		nnss <- strsplit(nss, split="\\.")[[1]]
		xpp <- seq_along(as.numeric(names(list.sim[[nss]]$mean)))
		xpp <- round(seq(xpp[1], xpp[length(xpp)], length.out=min(max.points, length(xpp))))
		lines(as.numeric(names(list.sim[[nss]]$mean))[xpp], sapply(list.sim[[nss]]$mean, function(x) FUN(x[[what]]))[xpp], col=col.sim[nnss[1]], lty=lty.sim[nnss[2]], pch=pch.sim[nnss[2]], type="o", lwd=lwd)
	}
}

list.sim <- list()

# Null model: no direct selection on robustness
list.sim$oo.o <- pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=10, grad.rob=c(0,0,0,0,0)), reps=reps, series.name="pure-null", force.run=force.run)

list.sim$ie.m <- pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=10, grad.rob=c(-grad.effect,0,0,0,0)), reps=reps, series.name="pure-ie-m", force.run=force.run)
list.sim$le.m <- pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=10, grad.rob=c(0,-grad.effect,0,0,0)), reps=reps, series.name="pure-le-m", force.run=force.run)
list.sim$im.m <- pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=10, grad.rob=c(0,0,-grad.effect,0,0)), reps=reps, series.name="pure-im-m", force.run=force.run)
list.sim$lm.m <- pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=10, grad.rob=c(0,0,0,-grad.effect,0)), reps=reps, series.name="pure-lm-m", force.run=force.run)
list.sim$st.m <- pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=10, grad.rob=c(0,0,0,0,-grad.effect)), reps=reps, series.name="pure-st-m", force.run=force.run)

list.sim$ie.p <- pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=10, grad.rob=c(grad.effect,0,0,0,0)), reps=reps, series.name="pure-ie-p", force.run=force.run)
list.sim$le.p <- pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=10, grad.rob=c(0,grad.effect,0,0,0)), reps=reps, series.name="pure-le-p", force.run=force.run)
list.sim$im.p <- pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=10, grad.rob=c(0,0,grad.effect,0,0)), reps=reps, series.name="pure-im-p", force.run=force.run)
list.sim$lm.p <- pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=10, grad.rob=c(0,0,0,grad.effect,0)), reps=reps, series.name="pure-lm-p", force.run=force.run)
list.sim$st.p <- pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=10, grad.rob=c(0,0,0,0,grad.effect)), reps=reps, series.name="pure-st-p", force.run=force.run)


pdf("figG.pdf", width=12, height=8)
layout(rbind(1:3,4:6))
par(cex=1)
allplots(list.sim, what="initenv")
allplots(list.sim, what="lateenv")
allplots(list.sim, what="initmut")
allplots(list.sim, what="latemut")
allplots(list.sim, what="stability")
dev.off()
