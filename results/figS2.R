#!/usr/bin/env Rscript

# Compute and plot correlations among robustness indexes 
# from random or evolved gene networks

source("./terminology.R")
source("./defaults.R")

source("../src/randnetwork.R")
source("../src/robindex.R")

################## Options

param <- default

param$rob.reps     <- 1000
net.reps           <- 10000

cache.tag <- "figA"

# Graphical options
maxplotpoints <- 1000 # avoids overcrowded plots
xylims        <- c(-40,-2)


################## Running simulations

cache.file <- paste0(param$cache.dir, "/", cache.tag, "-main.rds")

res <- if (param$use.cache && file.exists(cache.file)) {
			readRDS(cache.file)
		} else {
			Wlist <- replicate(net.reps, randW(net.size=param$n, reg.mean=param$rand.mean, reg.sd=param$rand.sd, density=param$density), simplify=FALSE)
			ans <- mclapply(Wlist, function(W) {
				robindex.Wmatrix(
				W               = W, 
				a               = param$a,
				dev.steps       = param$dev.steps,
				measure         = param$dev.measure, 
				mut.sd          = param$sim.mutsd, 
				mut.correlated  = param$mut.correlated, 
				test.initmut.sd = param$initenv.sd, 
				test.latemut.sd = param$lateenv.sd, 
				test.initenv.sd = param$initenv.sd,
				test.lateenv.sd = param$lateenv.sd, 
				test.reps       = param$rob.reps, 
				log.robustness  = param$log.robustness)
			}, mc.cores=param$mc.cores)
			saveRDS(ans, file=cache.file, version=2)
			ans
		}

##################### Figure

lp <- length(phen.expression)
mm <- matrix(0, ncol=lp-1, nrow=lp-1)
mm[lower.tri(mm, diag=TRUE)] <- 1:(lp*(lp-1)/2)

pdf(paste0("figS2.pdf"), width=param$maxfigwidth/param$figscale, height=param$maxfigwidth/param$figscale, pointsize=param$pointsize)
	layout(mm)
	par(mar=0.1+c(0,0,0,0), oma=c(4,5,0,0), cex=1)
	for (ii in 1:(lp-1)) {
	    for (jj in ((ii+1):lp)) {
	        rrx <- sapply(res, function(x) param$summary.FUN(x[[names(phen.expression)[ii]]]))
	        rry <- sapply(res, function(x) param$summary.FUN(x[[names(phen.expression)[jj]]]))
	        plot(rrx[1:min(length(res), maxplotpoints)], rry[1:min(length(res), maxplotpoints)], xaxt="n", yaxt="n", xlab="", ylab="", col="gray", xlim=xylims, ylim=xylims)
	        if (ii==1) {
	            axis(2)
	            mtext(as.expression(phen.expression[jj]), col=default.cols[jj], side=2, line=3)
	        }
	        if (jj==lp) {
	            axis(1)
	            mtext(as.expression(phen.expression[ii]), col=default.cols[ii], side=1, line=3)
	        }
	        legend("topleft", paste0("r=", format(round(cor(rrx, rry), digits=2), nsmall=2)), bty="n", cex=1)
	    }
	}
dev.off()

