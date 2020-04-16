#
# Evolvabilities and conditional evolvabilities (from M matrices)

source("./terminology.R")
source("./defaults.R")
source("./randnetwork.R")
source("../src/netw.R")
source("../src/robindex.R")

library(parallel)
mc.cores <- default.mc.cores

Wtoconsider <- c("random", "evolved", "randevol")
use.cache <- TRUE
cache.dir <- "../cache"

M.reps        <- 100
rob.reps      <- default.rob.reps
a             <- default.a
dev.steps     <- default.dev.steps
measure       <- default.dev.measure
initenv.sd    <- default.initenv.sd
lateenv.sd    <- default.lateenv.sd
initmut.sd    <- default.initmut.sd
latemut.sd    <- default.latemut.sd
mut.sd        <- default.sim.mutsd
mut.correlated <- default.mut.correlated

net.size      <- default.n
netsel.size   <- default.nsel

epsilon.zero  <- default.epsilon.zero

force.run     <- !use.cache

# For random matrices
W.reps        <- 100
rand.density  <- 1
rand.mean     <- -0.2
rand.sd       <- 1.2
min.evolv     <- 1e-3 # could be NA (no plotting)

# For evolved matrices
evolved.file.pattern <- 'figG-null-\\d+.rds'
evolved.gen          <- NA    # NA: last generation of the simulations

nb.mut  <- 1

boxplot.outline  <- FALSE

mutate <- function(W, mut.sd) {
	which.mut <- sample(size=1,  which(W != 0)) # Bug if only one W != 0
	W[which.mut] <- rnorm(1, mean=if(mut.correlated) W[which.mut] else 0, sd=mut.sd)
	W
}

fullPhen <- function(W) {
	phen <- model.M2 (W=W, a=a, steps=dev.steps)$mean
	initenv <- robindex.initenv(W=W, a=a, dev.steps=dev.steps, measure=measure, env.sd=initenv.sd, rep=rob.reps, log=TRUE)
	lateenv <- robindex.lateenv(W=W, a=a, dev.steps=dev.steps, measure=measure, env.sd=lateenv.sd, rep=rob.reps, log=TRUE)
	initmut <- robindex.initmut(W=W, a=a, dev.steps=dev.steps, measure=measure, mut.sd=initmut.sd, rep=rob.reps, log=TRUE)
	latemut <- robindex.latemut(W=W, a=a, dev.steps=dev.steps, measure=measure, mut.sd=latemut.sd, rep=rob.reps, log=TRUE)
	stability <- robindex.stability(W=W, a=a, dev.steps=dev.steps, measure=measure, log=TRUE)
	ans <- c(phen, initenv, lateenv, initmut, latemut, stability, mean(initenv), mean(lateenv), mean(initmut), mean(latemut), mean(stability))
	xx <- 1:nrow(W)
	names(ans) <- c(paste0("phen.", xx), paste0("initenv.", xx), paste0("lateenv.", xx), paste0("initmut.", xx), paste0("latemut.", xx), paste0("stability.", xx), "initenv.mean", "lateenv.mean", "initmut.mean", "latemut.mean", "stability.mean")
	ans
}

Mmatrices <- function(W, nbmut=nb.mut, reps=M.reps) {
	do.call(rbind, mclapply(1:reps, function(i) {
		myW <- W
		for (i in 1:nbmut)
			myW <- mutate(myW, mut.sd)
		fullPhen(myW)
		}, mc.cores=mc.cores)
	)
}

fullM <- function(W, nbmut=nb.mut) {
	var(Mmatrices(W, nbmut))
}

condEvolv <- function(G, focal=1, conditional=2:ncol(G)) {
	G <- G[c(focal, conditional),c(focal,conditional)]
	beta <- c(1, rep(0, length(conditional)))
	ans <- try(solve(t(beta)%*%solve(G)%*% beta), silent=TRUE)
	if (class(ans) == "try-error") min.evolv else ans
}

plotM <- function(W, trait1, trait2, nbmut=nb.mut, rep=1000, ...) {
	mm <- Mmatrices(W, nbmut=nbmut, reps=rep)
	plot(mm[,paste0(trait1, ".mean")], mm[,paste0(trait2, ".mean")], xlab=as.expression(phen.expression[trait1]), ylab=as.expression(phen.expression[trait2]), ...)
}

robs <- c("initenv", "lateenv", "initmut", "latemut", "stability")
names(robs) <- c(TERM.ENVCAN.SHORT, TERM.HOMEO.SHORT, TERM.GENCAN.SHORT, TERM.SOM.SHORT, TERM.STAB.SHORT)
cols <- c(initenv=COL.ENVCAN, lateenv=COL.HOMEO, initmut=COL.GENCAN, latemut=COL.SOM, stability=COL.STAB)


for (Wstyle in Wtoconsider) {

	cache.file <- paste0(cache.dir, "/figJ-", Wstyle, ".rds")
	if (!dir.exists(cache.dir)) dir.create(cache.dir)	

	dd <- NULL # list of M matrices
	dd <- if (use.cache && file.exists(cache.file)) readRDS(cache.file)
	
	evolved.files <- list.files(path=cache.dir, pattern=evolved.file.pattern, full.names=TRUE)
	
	reg.mean <- rand.mean
	reg.sd   <- rand.sd
	reg.density <- rand.density
	if (Wstyle == "randevol") {
		Wevoldist <- Wdist.fromfiles(evolved.files, epsilon.zero=epsilon.zero)
		reg.mean <- Wevoldist$mean
		reg.sd   <- Wevoldist$sd
		reg.density<- Wevoldist$density
	}	

	if (is.null(dd)) {
		dd <- mclapply(evolved.files, function(r) {
			if (Wstyle == "evolved") {
				ss <- readRDS(r)
				if (is.na(evolved.gen) || !as.character(evolved.gen) %in% names(ss)) evolved.gen <- names(ss)[length(ss)]
				W <- ss[[as.character(evolved.gen)]]$W
			} else { # both random and randevol
				W <- randW(net.size, reg.mean, reg.sd, reg.density)
			}
			return(fullM(W))
		}, mc.cores=mc.cores)
		saveRDS(dd, file=cache.file)
	}
	
	evolv.free  <- lapply(robs, function(rob) sapply(dd, function(M) condEvolv(M, paste0(rob, ".mean"), NULL)))
	evolv.phens <- lapply(robs, function(rob) sapply(dd, function(M) condEvolv(M, paste0(rob, ".mean"), paste0("phen.", 1:netsel.size))))
	evolv.robs  <- lapply(robs, function(rob) sapply(dd, function(M) condEvolv(M, paste0(rob, ".mean"), paste0(robs[!robs %in% rob], ".mean"))))
	
	if (sum(unlist(evolv.free) == min.evolv) > 0) warning(Wstyle, ":", sum(unlist(evolv.free) == min.evolv), "/", length(unlist(evolv.free)), " unconditional evolvabilities not computed (singular G)")
	if (sum(unlist(evolv.phens) == min.evolv) > 0) warning(Wstyle, ":", sum(unlist(evolv.phens) == min.evolv), "/", length(unlist(evolv.phens)), " phenotypic conditional evolvabilities not computed (singular G)")
	if (sum(unlist(evolv.robs) == min.evolv) > 0) warning(Wstyle, ":", sum(unlist(evolv.robs) == min.evolv), "/", length(unlist(evolv.robs)), " robustness conditional evolvabilities not computed (singular G)")

	
	lr <- length(robs)
	lt <- 3 # number of evolvabilities to display
	
	pdf(paste0("figJ-", Wstyle, ".pdf"), width=8, height=5)
		
		boxplot(evolv.free, at=1+(0:(lr-1))*(lt+1), xlim=c(0, (lt+1)*lr), log="y", xaxt="n",  ylab="Evolvability", border=cols, density=10, ylim=c(0.01,100), outline=boxplot.outline)
		boxplot(evolv.phens, at=2+(0:(lr-1))*(lt+1), xaxt="n", add=TRUE, border=cols, col="lightgray", outline=boxplot.outline)
		boxplot(evolv.robs, at=3+(0:(lr-1))*(lt+1), xaxt="n", add=TRUE, border=cols, col="bisque", outline=boxplot.outline)
		
		axis(1, at=1+(lt-1)/2+(0:(lr-1))*(lt+1), labels=names(robs), tick=FALSE)
		legend("topleft", pch=22, pt.bg=c("white","lightgray","bisque"), legend=c("Unconditional", "Cond. gene expression","Cond. other robustness"), horiz=TRUE, bty="n", pt.cex=1.5)
	
	dev.off()
}

#~ set.seed(2)
#~ example.Mrand <- randW.default(net.size=6, reg.mean=0, reg.sd=1)
#~ example.Wrand <- matrix(rnorm(36), ncol=6)
#~ example.Wsim <- readRDS("../cache/figG-null-1.rds")[['5000']]$W
#~ plotM(example.Wrand, trait1="latemut", trait2="stability", rep=1000, nbmut=3)
#~ plotM(example.Wsim, trait1="latemut", trait2="stability", rep=1000, nbmut=1)

