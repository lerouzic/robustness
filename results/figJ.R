#
# Evolvabilities and conditional evolvabilities (from M matrices)

source("terminology.R")

source("../src/netw.R")
source("../src/robindex.R")

library(parallel)
mc.cores <- min(128, detectCores()-1)

reps <- 10000
rob.reps <- 1000
a <- 0.2
dev.steps <- 16
initenv.sd <- 0.1
lateenv.sd <- 0.1
initmut.sd <- 0.1
latemut.sd <- 0.1
mut.sd     <- 0.1
mut.correlated <- TRUE
force.run <- FALSE

ref.sim <- "pure-null"
gen.sim <- "5000"

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

condEvolv <- function(G, focal=1, conditional=2:ncol(G)) {
	G <- G[c(focal, conditional),c(focal,conditional)]
	beta <- c(1, rep(0, length(conditional)))
	solve(t(beta)%*%solve(G)%*% beta)
}

robs <- c("initenv", "lateenv", "initmut", "latemut", "stability")
names(robs) <- c(TERM.ENVCAN.SHORT, TERM.HOMEO.SHORT, TERM.GENCAN.SHORT, TERM.SOM.SHORT, TERM.STAB.SHORT)
cols <- c(initenv=COL.ENVCAN, lateenv=COL.HOMEO, initmut=COL.GENCAN, latemut=COL.SOM, stability=COL.STAB)

Mmats.cachefile <- "../cache/Mmats.rds"

if (force.run || !file.exists(Mmats.cachefile)) {
	Wmats <- lapply(list.files(path="../cache", pattern=paste0(ref.sim, "-.*\\.rds"),full.names=TRUE), function(ff) readRDS(ff)[[gen.sim]]$W)
	Mmats <- lapply(Wmats, fullM)
	saveRds(Mmats, file=Mmats.cachefile)
}
Mmats <- readRDS(Mmats.cachefile)

evolv.free <- lapply(robs, function(rob) sapply(Mmats, function(M) condEvolv(M, paste0(rob, ".mean"), NULL)))
evolv.phen <- lapply(robs, function(rob) sapply(Mmats, function(M) condEvolv(M, paste0(rob, ".mean"), paste0("phen.", 1:ncol(Wmats[[1]])))))
evolv.phens <- lapply(robs, function(rob) sapply(Mmats, function(M) condEvolv(M, paste0(rob, ".mean"), paste0("phen.", 1:3))))
evolv.robs <- lapply(robs, function(rob) sapply(Mmats, function(M) condEvolv(M, paste0(rob, ".mean"), paste0(robs[!robs %in% rob], ".mean"))))

lr <- length(robs)
lt <- 3

pdf("figJ.pdf", width=8, height=5)

boxplot(evolv.free, at=1+(0:(lr-1))*(lt+1), xlim=c(0, (lt+1)*lr), log="y", xaxt="n",  ylab="Evolvability", border=cols, density=10, ylim=c(0.01,20))
boxplot(evolv.phens, at=2+(0:(lr-1))*(lt+1), xaxt="n", add=TRUE, border=cols, col="lightgray")
boxplot(evolv.robs, at=3+(0:(lr-1))*(lt+1), xaxt="n", add=TRUE, border=cols, col="bisque")

axis(1, at=1+(lt-1)/2+(0:(lr-1))*(lt+1), labels=names(robs), tick=FALSE)
legend("topleft", pch=22, bg=c("white","lightgray","bisque"), legend=c("Unconditional", "Cond. gene expression","Cond. other robustness"), horiz=TRUE, cex=1.5)

dev.off()
