#!/usr/bin/env Rscript

# Compute and plot correlations among robustness indexes 
# from random gene networks

source("./terminology.R")
source("../src/robindex.R")
library(parallel)
mc.cores <- min(detectCores()-1, 64)

phen <- c(
    #mean=TERM.EXPRESSION,
    initenv=substitute(x~(y), list(x=TERM.ENVCAN.LONG, y=ABBRV.ENVCAN[[1]])),
    lateenv=substitute(x~(y), list(x=TERM.HOMEO.LONG, y=ABBRV.HOMEO[[1]])),
    initmut=substitute(x~(y), list(x=TERM.GENCAN.LONG, y=ABBRV.GENCAN[[1]])),
    latemut=substitute(x~(y), list(x=TERM.SOM.LONG, y=ABBRV.SOM[[1]])),
    stability=substitute(x~(y), list(x=TERM.STAB.LONG, y=ABBRV.STAB[[1]])))

use.cache <- TRUE

net.size <- 10
reps <- 10000
density <- 0.2
reg.mean <- -0.2
reg.sd <- 1.2

a <- 0.2
dev.steps <- 20

rob.reps <- 100
rob.initenv.sd <- 0.1
rob.lateenv.sd  <- 0.1
rob.mut.sd     <- 0.1

maxplotpoints <- min(reps, 1000) # avoids overcrowded plots
lowthresh <- 1e-4 # this is necessary when plotting log variances

cache.dir <- "../cache"
cache.file <- paste(cache.dir, "figA.Rda", sep="/")
if (!dir.exists(cache.dir)) dir.create(cache.dir)

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
}

if (length(dd) < reps || nrow(dd[[1]]$W) != net.size)
	warning("The simulations stored in the cache do not fit the parameters. Consider cleaning the cache and run the script again.")

saveRDS(dd, cache.file)

lp <- length(phen)
mm <- matrix(0, ncol=lp-1, nrow=lp-1)
mm[lower.tri(mm, diag=TRUE)] <- 1:(lp*(lp-1)/2)

whattoconsider<- function(x) x[[1]] # the first gene of the network
whattoconsider <- function(x) mean(x) # the average index for all genes

pdf("figA.pdf", width=10, height=10)
layout(mm)
par(mar=0.1+c(0,0,0,0), oma=c(4,5,0,0))
for (ii in 1:(lp-1)) {
    for (jj in ((ii+1):lp)) {
        rrx <- sapply(dd[1:maxplotpoints], function(x) whattoconsider(x[[names(phen)[ii]]]))
        rry <- sapply(dd[1:maxplotpoints], function(x) whattoconsider(x[[names(phen)[jj]]]))
#~         rrx[rrx < lowthresh] <- lowthresh
#~         rry[rry < lowthresh] <- lowthresh
        plot(rrx, rry, xaxt="n", yaxt="n", xlab="", ylab="", col="gray")
        if (ii==1) {
            axis(2)
            mtext(as.expression(phen[jj]), side=2, line=3)
        }
        if (jj==lp) {
            axis(1)
            mtext(as.expression(phen[ii]), side=1, line=3)
        }
        # abline(lm( rr[,names(phen)[jj]] ~ rr[,names(phen)[ii]]), col="red")
        legend("topleft", paste0("r=", round(cor(rrx, rry), digits=2)), bty="n", cex=1.5)
    }
}
dev.off()

