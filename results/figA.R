#!/usr/bin/env Rscript

# Compute and plot correlations among robustness indexes 
# from random or evolved gene networks

source("./terminology.R")
source("./defaults.R")
source("./randnetwork.R")
source("../src/robindex.R")

library(parallel)
mc.cores <- default.mc.cores

use.cache <- TRUE
cache.dir <- "../cache"

phen <- phen.expression   # from terminology.R

Wstyle <- "random" # Possible: "random", "evolved", "randevol"
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
xylims        <- list(random=c(-40,-2), evolved=NULL, randevol=NULL) 

# For random matrices
reps           <- 10000
net.size       <- default.n
rand.density   <- default.density
rand.mean      <- default.rand.mean
rand.sd        <- default.rand.sd

# For evolved matrices
evolved.file.pattern <- 'figG-null-\\d+.rds'
evolved.gen          <- NA    # NA: last generation of the simulations

# For density estimates from evolved matrices
epsilon.zero    <- default.epsilon.zero # W values below this will be considered as zero

cache.file <- paste0(cache.dir, "/figA-", Wstyle, ".rds")
if (!dir.exists(cache.dir)) dir.create(cache.dir)	

dd <- NULL
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
	dd <- mclapply(if (Wstyle=="evolved") evolved.files else seq_len(reps), function(r) {
		if (Wstyle == "evolved") {
			ss <- readRDS(r)
			if (is.na(evolved.gen) || !as.character(evolved.gen) %in% names(ss)) evolved.gen <- names(ss)[length(ss)]
			W <- ss[[as.character(evolved.gen)]]$W
		} else { # both random and randevol
			W <- randW(net.size, reg.mean, reg.sd, reg.density)
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
	saveRDS(dd, cache.file)
}

lp <- length(phen)
mm <- matrix(0, ncol=lp-1, nrow=lp-1)
mm[lower.tri(mm, diag=TRUE)] <- 1:(lp*(lp-1)/2)

pdf(paste0("figA.pdf"), width=10, height=10)
	layout(mm)
	par(mar=0.1+c(0,0,0,0), oma=c(4,5,0,0))
	for (ii in 1:(lp-1)) {
	    for (jj in ((ii+1):lp)) {
	        rrx <- sapply(dd[1:min(length(dd), maxplotpoints)], function(x) whattoconsider(x[[names(phen)[ii]]]))
	        rry <- sapply(dd[1:min(length(dd), maxplotpoints)], function(x) whattoconsider(x[[names(phen)[jj]]]))
	        plot(rrx, rry, xaxt="n", yaxt="n", xlab="", ylab="", col="gray", xlim=xylims[[Wstyle]], ylim=xylims[[Wstyle]])
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

