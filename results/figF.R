#!/usr/bin/env Rscript

# Evolution of robustness G matrices

source("./commonsim.R")
source("../src/robindex.R")

tags <- c("initenv", "lateenv", "initmut", "latemut", "stability")

library(parallel)
mc.cores <- min(64, detectCores()-1)

force.run=FALSE

#~ reps <- 24
#~ G <- 10000
#~ N <- 5000
#~ summary.every <- 1000

###### Specific for tests
reps <- 2
G <- 100
N <- 1000
summary.every <- 10
force.run <- TRUE
######

nselgen <- 3
new.s <- paste(c(rep(default.args$s, nselgen), rep(0, default.args$n-nselgen)), collapse=" ")
randtheta <- function() paste(round(runif(default.args$n, 0, 1), digits=3), collapse=" ")
stab <- sim.run.reps(args <- list(G=G, N=N, summary.every=summary.every, theta=randtheta(), s=new.s),reps=reps, series.name="stabtest", force.run=force.run, mc.cores=mc.cores)
#~ istabty <- sim.run.reps(args <- list(G=G, N=N, summary.every=summary.every, theta=randtheta(), s=new.s, ss=36000), reps=reps, series.name="stabty", force.run=force.run, mc.cores=mc.cores)
#~ som <- sim.run.reps(args <- list(G=G, N=N, summary.every=summary.every, theta=randtheta(), s=new.s, som.rate=1, som.sd=0.1),reps=reps, series.name="som", force.run=force.run, mc.cores=mc.cores)
#~ ienv <- sim.run.reps(args <- list(G=G, N=N, summary.every=summary.every, theta=randtheta(), s=new.s, initenv.sd=0.1), reps=reps, series.name="ienv", force.run=force.run, mc.cores=mc.cores)
#~ lenv <- sim.run.reps(args <- list(G=G, N=N, summary.every=summary.every, theta=randtheta(), s=new.s, lateenv.sd=0.1), reps=reps, series.name="lenv", force.run=force.run, mc.cores=mc.cores)


# Individual sims
for (ssi in seq_along(stab$full))
	plotmat(robindex.Gmatrix.outfile(stab$full[[ssi]]), arg1=2, arg2=3, col="blue", col.end="green", add=ssi!=1)

for (ssi in seq_along(stab$full))
	plotmat(do.call(robindex.Mmatrix.outfile, c(list(out=stab$full[[ssi]], rep=100), default.args[names(default.args) %in% names(formals(robindex.Mmatrix))])), 
		arg1=2, arg2=3, col="blue", col.end="green", add=ssi!=1)


# All sims
plotmat(robindex.Gmatrix.outfile(stab$mean), arg1=2, arg2=3, col="blue", col.end="green", var.thresh=1e-20)

