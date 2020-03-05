#!/usr/bin/env Rscript

# Ensures that the last PCs are real (not due to stochastic noise)

source("./terminology.R")
source("../src/robindex.R")
library(parallel)
mc.cores <- min(detectCores()-1, 64)

cache.dir <- "../cache"
if (!dir.exists(cache.dir)) dir.create(cache.dir)

use.cache <- TRUE

cols <- 1:5

a <- 0.2
dev.steps <- 16
rob.initenv.sd <- 0.1
rob.lateenv.sd  <- 0.1
rob.mut.sd     <- 0.1

default.reg.mean <- -0.2
default.reg.sd   <-  1.2
default.net.size <- 10
default.density  <- 1
default.reps <- 1000
default.rob.reps <- 100

#~ whattoconsider<- function(x) x[[1]] # the first gene of the network
whattoconsider <- function(x) mean(x) # the average index for all genes

eigenV <- function(reps, rob.reps, reg.mean=default.reg.mean, reg.sd=default.reg.sd, net.size=default.net.size, density=default.density) {
	
	cache.file <- paste0(cache.dir, "/", "figH-R", reps, "-r", rob.reps, ".Rda")
	dd <- NULL
	dd <- if (use.cache && file.exists(cache.file)) readRDS(cache.file)
	if (is.null(dd)) {
		dd <- mclapply(1:reps, function(r) {
			W <- matrix(rnorm(net.size^2, mean=reg.mean, sd=reg.sd), ncol=net.size)
			W[sample.int(net.size^2, floor((1-density)*net.size^2))] <- 0
			list(W=W, 
				mean=model.M2(W, a, steps=dev.steps)$mean, 
				initenv=robindex.initenv(W, a, dev.steps, rob.initenv.sd, rep=rob.reps, log=TRUE),
				lateenv=robindex.lateenv(W, a, dev.steps, rob.lateenv.sd, rep=rob.reps, log=TRUE),
				initmut=robindex.initmut(W, a, dev.steps, rob.mut.sd, rep=rob.reps,log=TRUE),
				latemut=robindex.latemut(W, a, dev.steps, rob.mut.sd, rep=rob.reps, log=TRUE),
				stability=robindex.stability(W, a, dev.steps, log=TRUE)
			)
		}, mc.cores=mc.cores) 
		saveRDS(dd, cache.file)
	}
	
	rrr <- do.call(rbind, lapply(dd, function(ddd) sapply(c("initenv", "lateenv", "initmut", "latemut", "stability"), function(ppp) whattoconsider(ddd[[ppp]]))))
	prp <- prcomp(rrr, scale.=TRUE)
	prp$sdev^2/(sum(prp$sdev^2))
}

allreps <- round(10^(seq(2, 4, length.out=7)))
resreps <- lapply(allreps, function(reps) eigenV(reps, default.rob.reps))

allrobs <- round(10^(seq(1, 4, length.out=7)))
resrobs <- lapply(allrobs, function(robs) eigenV(default.reps, robs))

pdf("figH.pdf", width=5, height=10)
layout(1:2)

plot(NULL, xlim=range(allreps), ylim=c(1e-2,1), log="xy", xlab="Number of simulated networks", ylab="Proportion variance explained by each PC")
for (i in 1:5)
	lines(allreps, sapply(resreps, function(r) r[i]), col=cols[i])
text(allreps[1], resreps[[1]], col=cols, pos=1, paste0("PC", seq_along(cols)))

plot(NULL, xlim=range(allrobs), ylim=c(1e-2,1), log="xy", xlab="Number of robustness tests", ylab="Proportion variance explained by each PC")
for (i in 1:5)
	lines(allrobs, sapply(resrobs, function(r) r[i]), col=cols[i])
text(allrobs[1], resrobs[[1]], col=cols, pos=1, paste0("PC", seq_along(cols)))	
	
dev.off()
