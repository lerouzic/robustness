source("commonpure.R")
source("terminology.R")

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
lty.sim <- c(p=2, m=1, o=3)

allplots <- function(list.sim, what="fitness", xlim=NULL, ylim=NULL, xlab="Generations", ylab=what, lwd=3, FUN=mean, ...) {
	if(is.null(xlim)) xlim <- c(0, as.numeric(rev(names(list.sim[[1]]$mean))[1]))
	allv <- do.call(rbind, lapply(list.sim, function(x) do.call(rbind, lapply(x$mean, function(xx) FUN(xx[[what]])))))
	if(is.null(ylim)) ylim <- c(min(allv), max(allv))
	plot(NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
	for (nss in names(list.sim)) {
		nnss <- strsplit(nss, split="\\.")[[1]]
		lines(as.numeric(names(list.sim[[nss]]$mean)), sapply(list.sim[[nss]]$mean, function(x) FUN(x[[what]])), col=col.sim[nnss[1]], lty=lty.sim[nnss[2]], lwd=lwd)
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


allplots(list.sim, what="initenv")
allplots(list.sim, what="lateenv")
allplots(list.sim, what="initmut")
allplots(list.sim, what="latemut")
allplots(list.sim, what="stability")
