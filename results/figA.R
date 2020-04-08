#!/usr/bin/env Rscript

# Compute and plot correlations among robustness indexes 
# from random or evolved gene networks

source("./terminology.R")
source("./defaults.R")
source("../src/robindex.R")

library(parallel)
mc.cores <- default.mc.cores
mc.cores <- 1

use.cache <- TRUE
cache.dir <- "../cache"

phen <- phen.expression   # from terminology.R

Wtoconsider <- c("random", "evolved")
# whattoconsider<- function(x) x[[1]] # the first gene of the network
whattoconsider <- function(x) mean(x) # the average index for all genes

# Most parameters use the defaults, some specificities
a              <- default.a
dev.steps      <- default.dev.steps
measure        <- default.dev.measure
log.robustness <- default.log.robustness

rob.reps       <- 1000
rob.initenv.sd <- default.initenv.sd
rob.lateenv.sd <- default.lateenv.sd
rob.initmut.sd <- default.initmut.sd
rob.latemut.sd <- default.latemut.sd

# Graphical options
maxplotpoints <- 1000 # avoids overcrowded plots
xylims        <- c(-40,-2) # can be NULL

# For random matrices
net.size       <- 6
reps           <- 10000
density        <- 1
reg.mean       <- -0.2
reg.sd         <- 1.2

# For evolved matrices
evolved.file.pattern <- 'figG-null-\\d+.rds'
evolved.gen          <- NA    # NA: last generation of the simulations



for (Wstyle in Wtoconsider) {

	cache.file <- paste0(cache.dir, "/figA-", Wstyle, ".rds")
	if (!dir.exists(cache.dir)) dir.create(cache.dir)	

	dd <- NULL
	dd <- if (use.cache && file.exists(cache.file)) readRDS(cache.file)
	
	loopover <- if (Wstyle == "random") 1:reps else if (Wstyle == "evolved") list.files(path=cache.dir, pattern=evolved.file.pattern, full.names=TRUE)

	if (is.null(dd)) {
		dd <- mclapply(loopover, function(r) {
			if (Wstyle == "random") {
				W <- matrix(rnorm(net.size^2, mean=reg.mean, sd=reg.sd), ncol=net.size)
				W[sample.int(net.size^2, floor((1-density)*net.size^2))] <- 0
			} else if (Wstyle == "evolved") {
				ss <- readRDS(r)
				if (is.na(evolved.gen) || !as.character(evolved.gen) %in% names(ss)) evolved.gen <- names(ss)[length(ss)]
				W <- ss[[as.character(evolved.gen)]]$W
			}
			
			list(W=W, 
				mean=model.M2(W, a, steps=dev.steps, measure=measure)$mean, 
				initenv=robindex.initenv(W, a, dev.steps, measure, rob.initenv.sd, rep=rob.reps, log=log.robustness),
				lateenv=robindex.lateenv(W, a, dev.steps, measure, rob.lateenv.sd, rep=rob.reps, log=log.robustness),
				initmut=robindex.initmut(W, a, dev.steps, measure, rob.initmut.sd, rep=rob.reps, log=log.robustness),
				latemut=robindex.latemut(W, a, dev.steps, measure, rob.latemut.sd, rep=rob.reps, log=log.robustness),
				stability=robindex.stability(W, a, dev.steps, measure, log=log.robustness)
			)
		}, mc.cores=mc.cores) 
	}

	saveRDS(dd, cache.file)

	lp <- length(phen)
	mm <- matrix(0, ncol=lp-1, nrow=lp-1)
	mm[lower.tri(mm, diag=TRUE)] <- 1:(lp*(lp-1)/2)

	pdf(paste0("figA-", Wstyle, ".pdf"), width=10, height=10)
		layout(mm)
		par(mar=0.1+c(0,0,0,0), oma=c(4,5,0,0))
		for (ii in 1:(lp-1)) {
		    for (jj in ((ii+1):lp)) {
		        rrx <- sapply(dd[1:min(length(dd), maxplotpoints)], function(x) whattoconsider(x[[names(phen)[ii]]]))
		        rry <- sapply(dd[1:min(length(dd), maxplotpoints)], function(x) whattoconsider(x[[names(phen)[jj]]]))
		        plot(rrx, rry, xaxt="n", yaxt="n", xlab="", ylab="", col="gray", xlim=xylims, ylim=xylims)
		        if (ii==1) {
		            axis(2)
		            mtext(as.expression(phen[jj]), side=2, line=3)
		        }
		        if (jj==lp) {
		            axis(1)
		            mtext(as.expression(phen[ii]), side=1, line=3)
		        }
		        # abline(lm( rr[,names(phen)[jj]] ~ rr[,names(phen)[ii]]), col="red")
		        legend("topleft", paste0("r=", format(round(cor(rrx, rry), digits=2), nsmall=2)), bty="n", cex=1.5)
		    }
		}
	dev.off()
 } # end of the for loop over Wstyles
