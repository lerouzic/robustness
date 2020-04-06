#!/usr/bin/env Rscript

# Ensures that the last PCs are real (not due to stochastic noise)

source("./terminology.R")
source("./defaults.R")
source("../src/robindex.R")

library(parallel)
mc.cores <- default.mc.cores

cache.dir <- "../cache"
if (!dir.exists(cache.dir)) dir.create(cache.dir)

use.cache <- TRUE

a               <- default.a        
dev.steps       <- default.dev.steps
measure         <- default.dev.measure
rob.initenv.sd  <- default.initenv.sd
rob.lateenv.sd  <- default.lateenv.sd
rob.initmut.sd  <- default.initmut.sd
rob.latemut.sd  <- default.latemut.sd  

# These should match figA and figB
reg.mean <- -0.2
reg.sd   <-  1.2
net.size <- default.n
density  <- 1

reps <- 5000
rob.reps <- 500

#~ whattoconsider<- function(x) x[[1]] # the first gene of the network
whattoconsider <- function(x) mean(x) # the average index for all genes

cols <- 1:5

eigenV <- function(reps, rob.reps) {
	# needs global variables: reg.mean, reg.sd, net.size, density
	cache.file <- paste0(cache.dir, "/", "figH-R", reps, "-r", rob.reps, ".rds")
	dd <- NULL
	dd <- if (use.cache && file.exists(cache.file)) readRDS(cache.file)
	if (is.null(dd)) {
		dd <- mclapply(1:reps, function(r) {
			W <- matrix(rnorm(net.size^2, mean=reg.mean, sd=reg.sd), ncol=net.size)
			W[sample.int(net.size^2, floor((1-density)*net.size^2))] <- 0
			list(W=W, 
				mean=model.M2(W, a, steps=dev.steps)$mean, 
				initenv=robindex.initenv(W, a, dev.steps, measure=measure, env.sd=rob.initenv.sd, rep=rob.reps, log=TRUE),
				lateenv=robindex.lateenv(W, a, dev.steps, measure=measure, env.sd=rob.lateenv.sd, rep=rob.reps, log=TRUE),
				initmut=robindex.initmut(W, a, dev.steps, measure=measure, mut.sd=rob.mut.sd, rep=rob.reps,log=TRUE),
				latemut=robindex.latemut(W, a, dev.steps, measure=measure, mut.sd=rob.mut.sd, rep=rob.reps, log=TRUE),
				stability=robindex.stability(W, a, dev.steps, measure=measure, log=TRUE)
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
allreps <- round(10^(seq(2, 4, length.out=7)))
resreps <- lapply(allreps, function(reps) eigenV(reps, default.rob.reps))

allrobs <- round(10^(seq(1, 5, length.out=9)))
resrobs <- lapply(allrobs, function(robs) eigenV(default.reps, robs))

pdf("figH.pdf", width=8, height=4)
	layout(t(1:2))
	
	plot(NULL, xlim=c(0.4,1)*range(allreps), ylim=c(1e-2,1), log="xy", xlab="Number of simulated networks", ylab="Proportion variance explained")
	for (i in 1:5)
		lines(allreps, sapply(resreps, function(r) r[i]), col=cols[i], type="o", pch=16)
	text(allreps[1], resreps[[1]], col=cols, pos=2, paste0("PC", seq_along(cols)))
	
	plot(NULL, xlim=c(0.2,1)*range(allrobs), ylim=c(1e-2,1), log="xy", xlab="Number of robustness tests", ylab="Proportion variance explained")
	for (i in 1:5)
		lines(allrobs, sapply(resrobs, function(r) r[i]), col=cols[i], type="o", pch=16)
	text(allrobs[1], resrobs[[1]], col=cols, pos=2, paste0("PC", seq_along(cols)))	
		
dev.off()
