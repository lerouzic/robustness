#
# Evolvabilities and conditional evolvabilities (from M matrices)

source("../src/netw.R")
source("../src/robindex.R")

library(parallel)
mc.cores <- min(128, detectCores()-1)

reps <- 1000
rob.reps <- 100
a <- 0.2
dev.steps <- 16
initenv.sd <- 0.1
lateenv.sd <- 0.1
initmut.sd <- 0.1
latemut.sd <- 0.1
mut.sd     <- 0.1
mut.correlated <- TRUE

mutate <- function(W, mut.sd) {
	which.mut <- sample(size=1,  which(W != 0)) # Bug if only one W != 0
	W[which.mut] <- rnorm(1, mean=if(mut.correlated) W[which.mut] else 0, sd=mut.sd)
	W
}

fullPhen <- function(W) {
	phen <- model.M2 (W=W, a=a, steps=dev.steps)$mean
	initenv <- robindex.initenv(W=W, a=a, dev.steps=dev.steps, env.sd=initenv.sd, rep=rob.reps, log=TRUE)
	lateenv <- robindex.lateenv(W=W, a=a, dev.steps=dev.steps, env.sd=lateenv.sd, rep=rob.reps, log=TRUE)
	initmut <- robindex.initmut(W=W, a=a, dev.steps=dev.steps, mut.sd=initmut.sd, rep=rob.reps, log=TRUE)
	latemut <- robindex.latemut(W=W, a=a, dev.steps=dev.steps, mut.sd=latemut.sd, rep=rob.reps, log=TRUE)
	stability <- robindex.stability(W=W, a=a, dev.steps=dev.steps, log=TRUE)
	ans <- c(phen, initenv, lateenv, initmut, latemut, stability, mean(initenv), mean(lateenv), mean(initmut), mean(latemut), mean(stability))
	xx <- 1:nrow(W)
	names(ans) <- c(paste0("phen.", xx), paste0("initenv.", xx), paste0("lateenv.", xx), paste0("initmut.", xx), paste0("latemut.", xx), paste0("stability.", xx), "initenv.mean", "lateenv.mean", "initmut.mean", "latemut.mean", "stability.mean")
	ans
}

fullM <- function(W) {
	mm <- do.call(rbind, mclapply(1:reps, function(i) {
		myW <- mutate(W, mut.sd)
		fullPhen(myW)
		}, mc.cores=mc.cores))
	var(mm)
}

W1 <- readRDS("../cache/pure-null-1.rds")[["5000"]]$W
W1.M <- fullM(W1)
