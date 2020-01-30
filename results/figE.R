#!/usr/bin/env Rscript

source("./commonsim.R")

library(parallel)
mc.cores <- min(64, detectCores()-1)

force.run=FALSE

reps <- 24
G <- 10000
N <- 5000
summary.every <- 1000
nselgen <- 3
new.s <- paste(c(rep(default.args$s, nselgen), rep(0, default.args$n-nselgen)), collapse=" ")
randtheta <- function() paste(round(runif(default.args$n, 0, 1), digits=3), collapse=" ")
stab <- sim.run.reps(args <- list(G=G, N=N, summary.every=summary.every, theta=randtheta(), s=new.s),reps=reps, series.name="stab", force.run=force.run, mc.cores=mc.cores)
istabty <- sim.run.reps(args <- list(G=G, N=N, summary.every=summary.every, theta=randtheta(), s=new.s, ss=36000), reps=reps, series.name="stabty", force.run=force.run, mc.cores=mc.cores)
som <- sim.run.reps(args <- list(G=G, N=N, summary.every=summary.every, theta=randtheta(), s=new.s, som.rate=1, som.sd=0.1),reps=reps, series.name="som", force.run=force.run, mc.cores=mc.cores)
ienv <- sim.run.reps(args <- list(G=G, N=N, summary.every=summary.every, theta=randtheta(), s=new.s, initenv.sd=0.1), reps=reps, series.name="ienv", force.run=force.run, mc.cores=mc.cores)
lenv <- sim.run.reps(args <- list(G=G, N=N, summary.every=summary.every, theta=randtheta(), s=new.s, lateenv.sd=0.1), reps=reps, series.name="lenv", force.run=force.run, mc.cores=mc.cores)

plotts(stab$mean, stab$sd, colname="fitness.mean", lwd=3, ylim=c(0,1))
plotts(stabty$mean, stabty$sd, colname="fitness.mean", lwd=3, add=TRUE, col="blue")
plotts(som$mean, som$sd, colname="fitness.mean", lwd=3, add=TRUE, col="brown")
plotts(ienv$mean, ienv$sd, colname="fitness.mean", lwd=3, add=TRUE, col="orange")
plotts(lenv$mean, lenv$sd, colname="fitness.mean", lwd=3, add=TRUE, col="darkgreen")

plotts(stab$mean, stab$sd, colname="robustness.initmut.mean", lwd=3, ylim=c(0,0.004))
plotts(stabty$mean, stabty$sd, colname="robustness.initmut.mean", lwd=3, add=TRUE, col="blue")
plotts(som$mean, som$sd, colname="robustness.initmut.mean", lwd=3, add=TRUE, col="brown")
plotts(ienv$mean, ienv$sd, colname="robustness.initmut.mean", lwd=3, add=TRUE, col="orange")
plotts(lenv$mean, lenv$sd, colname="robustness.initmut.mean", lwd=3, add=TRUE, col="darkgreen")

plotts(stab$mean, stab$sd, colname="robustness.latemut.mean", lwd=3, ylim=c(0,0.0015))
plotts(stabty$mean, stabty$sd, colname="robustness.latemut.mean", lwd=3, add=TRUE, col="blue")
plotts(som$mean, som$sd, colname="robustness.latemut.mean", lwd=3, add=TRUE, col="brown")
plotts(ienv$mean, ienv$sd, colname="robustness.latemut.mean", lwd=3, add=TRUE, col="orange")
plotts(lenv$mean, lenv$sd, colname="robustness.latemut.mean", lwd=3, add=TRUE, col="darkgreen")

plotts(stab$mean, stab$sd, colname="robustness.initenv.mean", lwd=3, ylim=c(0,0.0001))
plotts(stabty$mean, stabty$sd, colname="robustness.initenv.mean", lwd=3, add=TRUE, col="blue")
plotts(som$mean, som$sd, colname="robustness.initenv.mean", lwd=3, add=TRUE, col="brown")
plotts(ienv$mean, ienv$sd, colname="robustness.initenv.mean", lwd=3, add=TRUE, col="orange")
plotts(lenv$mean, lenv$sd, colname="robustness.initenv.mean", lwd=3, add=TRUE, col="darkgreen")

plotts(stab$mean, stab$sd, colname="robustness.lateenv.mean", lwd=3, ylim=c(0,0.003))
plotts(stabty$mean, stabty$sd, colname="robustness.lateenv.mean", lwd=3, add=TRUE, col="blue")
plotts(som$mean, som$sd, colname="robustness.lateenv.mean", lwd=3, add=TRUE, col="brown")
plotts(ienv$mean, ienv$sd, colname="robustness.lateenv.mean", lwd=3, add=TRUE, col="orange")
plotts(lenv$mean, lenv$sd, colname="robustness.lateenv.mean", lwd=3, add=TRUE, col="darkgreen")

plotts(stab$mean, stab$sd, colname="robustness.stability.mean", lwd=3, ylim=c(0,0.0001))
plotts(stabty$mean, stabty$sd, colname="robustness.stability.mean", lwd=3, add=TRUE, col="blue")
plotts(som$mean, som$sd, colname="robustness.stability.mean", lwd=3, add=TRUE, col="brown")
plotts(ienv$mean, ienv$sd, colname="robustness.stability.mean", lwd=3, add=TRUE, col="orange")
plotts(lenv$mean, lenv$sd, colname="robustness.stability.mean", lwd=3, add=TRUE, col="darkgreen")

