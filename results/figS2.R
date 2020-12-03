#!/usr/bin/env Rscript

# Ensures that the last PCs are real (not due to stochastic noise)

source("./terminology.R")
source("./defaults.R")
source("./randnetwork.R")
source("../src/robindex.R")

library(parallel)
mc.cores <- default.mc.cores

cache.dir <- "../cache"
if (!dir.exists(cache.dir)) dir.create(cache.dir)

use.cache <- TRUE
Wstyle <- "random" # Possible: "random", "evolved", "randevol"

a               <- default.a        
dev.steps       <- default.dev.steps
measure         <- default.dev.measure
rob.initenv.sd  <- default.initenv.sd
rob.lateenv.sd  <- default.lateenv.sd
rob.initmut.sd  <- default.initmut.sd
rob.latemut.sd  <- default.latemut.sd  
log.robustness  <- default.log.robustness

net.size        <- default.n
epsilon.zero    <- default.epsilon.zero

# These should match figA and figB
rand.mean       <- default.rand.mean
rand.sd         <- default.rand.sd
rand.density    <- default.density

reps            <- 5000
rob.reps        <- 500

allreps         <- round(10^(seq(2, 4, length.out=7)))
allrobs         <- round(10^(seq(1, 5, length.out=9)))

evolved.file.pattern <- 'figG-null-\\d+.rds'
evolved.files <- list.files(path=cache.dir, pattern=evolved.file.pattern, full.names=TRUE)
evolved.gen          <- NA

#~ whattoconsider<- function(x) x[[1]] # the first gene of the network
whattoconsider <- function(x) mean(x) # the average index for all genes

cols <- 1:5

eigenV <- function(reps, rob.reps, Wstyle) {
	# needs global variables: reg.mean, reg.sd, net.size, density
	cache.file <- paste0(cache.dir, "/", "figH-", Wstyle, "-R", reps, "-r", rob.reps, ".rds")
	
	dd <- NULL
	dd <- if (use.cache && file.exists(cache.file)) readRDS(cache.file)
	
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
		dd <- mclapply(evolved.files, function(f) {
			if (Wstyle == "evolved") {
				ss <- readRDS(f)
				if (is.na(evolved.gen) || !as.character(evolved.gen) %in% names(ss)) evolved.gen <- names(ss)[length(ss)]
				W <- ss[[as.character(evolved.gen)]]$W				
			} else {
				W <- randW(net.size, reg.mean, reg.sd, reg.density)
			}
			list(W=W, 
				mean=model.M2(W, a, steps=dev.steps)$mean, 
				initenv=robindex.initenv(W, a, dev.steps, measure=measure, env.sd=rob.initenv.sd, rep=rob.reps, log=log.robustness),
				lateenv=robindex.lateenv(W, a, dev.steps, measure=measure, env.sd=rob.lateenv.sd, rep=rob.reps, log=log.robustness),
				initmut=robindex.initmut(W, a, dev.steps, measure=measure, mut.sd=rob.initmut.sd, rep=rob.reps,log=log.robustness),
				latemut=robindex.latemut(W, a, dev.steps, measure=measure, mut.sd=rob.latemut.sd, rep=rob.reps, log=log.robustness),
				stability=robindex.stability(W, a, dev.steps, measure=measure, log=log.robustness)
			)
		}, mc.cores=mc.cores) 
		saveRDS(dd, cache.file)
	}
	rrr <- do.call(rbind, lapply(dd, function(ddd) sapply(c("initenv", "lateenv", "initmut", "latemut", "stability"), function(ppp) whattoconsider(ddd[[ppp]]))))
	prp <- prcomp(rrr, scale.=TRUE)
	prp$sdev^2/(sum(prp$sdev^2))
}

# Technically, this sounds pretty useless: it would be way more efficient to run the max number of replicates and sample them for lower
# counts. Yet, the current code is simpler, and avoids correlations among samples. 
resreps <- lapply(allreps, function(rreps) eigenV(rreps, default.rob.reps, Wstyle))
resrobs <- lapply(allrobs, function(robs) eigenV(reps, robs, Wstyle))

pdf(paste0("figS2.pdf"), width=8, height=4)
	layout(t(1:2))
	
	plot(NULL, xlim=c(0.4,1)*range(allreps), ylim=c(1e-3,1), log="xy", xlab="Number of simulated networks", ylab="Proportion variance explained")
	for (i in 1:5)
		lines(allreps, sapply(resreps, function(r) r[i]), col=cols[i], type="o", pch=16)
	text(allreps[1], resreps[[1]], col=cols, pos=2, paste0("PC", seq_along(cols)))
	
	plot(NULL, xlim=c(0.2,1)*range(allrobs), ylim=c(1e-3,1), log="xy", xlab="Number of robustness tests", ylab="Proportion variance explained")
	for (i in 1:5)
		lines(allrobs, sapply(resrobs, function(r) r[i]), col=cols[i], type="o", pch=16)
	text(allrobs[1], resrobs[[1]], col=cols, pos=2, paste0("PC", seq_along(cols)))	
	
dev.off()

