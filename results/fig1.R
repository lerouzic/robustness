#!/usr/bin/env Rscript

# Runs a PCA among robustness indexes and plots the results

source("./terminology.R")
source("./defaults.R")
source("../src/robindex.R") 
source("../src/randnetwork.R")
source("../src/tools.R")


#################### Options

Wstyle <- "random" # Possible: "random", "evolved" ,"randevol"
whattoconsider <- function(x) mean(x) # the average index for all genes

phen <- phen.expression      # from terminology.R    
phen.pos.ref <- "initmut"    # this guy will always be positive (so that the figure is reproducible)

################### Functions

plot.PCarrows <- function(pr, labels=default.labels, cols=default.cols, ...) {
	nPC <- length(pr$sdev)
	# In prcomp, the sign of the PCs is arbitrary. To ensure reproducibility, PCs will be oriented in such a was that phen.pos.ref > 0
    pr$rotation <- t(t(pr$rotation)*sign(pr$rotation[phen.pos.ref,]))
    
    plot(NULL, xlab="", ylab="", ylim=c(0.75,nPC+0.25), xlim=range(pr$rotation), yaxt="n", bty="n", xpd=NA, ...)
    arrows(x0=min(pr$rotation), y0=1:nPC, x1=max(pr$rotation), code=3, length=0.2)
    for (i in 1:nPC) text(x=pr$rotation[,i], y=nPC-i+1+0.1*(0:(nPC-1)), labels[rownames(pr$rotation)], col=cols[rownames(pr$rotation)], pos=3, cex=0.9, font=2, xpd=NA)
    lines(x=rep(0,2), y=c(1,nPC), lty=2)
    axis(4, at=1:nPC, labels=paste0("PC", nPC:1), las=2, tick=FALSE)
}

plot.PCbarplot <- function(pr, ...) {
	nPC <- length(pr$sdev)
	
	bb <- barplot(100*rev(pr$sdev^2/(sum(pr$sdev^2))), horiz=TRUE, ylim=c(1, 2*nPC), ylab="", xlab="% variance", width=1, space=1, ...)
    arrows(x0=seq(0, 80, 20), y0=1, y1=max(bb), col="gray", lty=3, length=0)
}

percentPC <- function(reps=1000, net.size=default.n, density=default.density, reg.mean=default.rand.mean, reg.sd=default.rand.sd, PC=1) {
	cache.str <- paste("figA", net.size, density, reg.mean, reg.sd, sep="-")
	cache.file <- paste0(cache.dir, "/", cache.str, ".rds")
	all.Wrob <-  
		if (file.exists(cache.file)) readRDS(cache.file)
		else {
			allW <- lapply(1:reps, function(i) randW(net.size=net.size, reg.mean=reg.mean, reg.sd=reg.sd, density=density))
			ans <- mclapply(allW, function(W) {
				list(
					initenv=robindex.initenv(W, 
						default.a, default.dev.steps, default.dev.measure, default.initenv.sd, rep=default.rob.reps, log=default.log.robustness),
					lateenv=robindex.lateenv(W, 
						default.a, default.dev.steps, default.dev.measure, default.lateenv.sd, rep=default.rob.reps, log=default.log.robustness),
					initmut=robindex.initmut(W, 
						default.a, default.dev.steps, default.dev.measure, default.initmut.sd, rep=default.rob.reps, log=default.log.robustness),
					latemut=robindex.latemut(W, 
						default.a, default.dev.steps, default.dev.measure, default.latemut.sd, rep=default.rob.reps, log=default.log.robustness),
					stability=robindex.stability(W, 
						default.a, default.dev.steps, default.dev.measure, log=default.log.robustness))
				}, mc.cores=mc.cores)
			saveRDS(ans, file=cache.file, version=2)
			ans
		}
	rrr <- do.call(rbind, lapply(all.Wrob, function(ddd) sapply(names(phen), function(ppp) whattoconsider(ddd[[ppp]]))))
	prp <- try(prcomp(rrr[, apply(rrr, 2, var) != 0], scale.=TRUE))
	if (class(prp) == "try-error") NA else (prp$sdev^2 / sum(prp$sdev^2))[PC]
}

#################### Calc
source("./figS2.R") # This script uses dataset figA. Run figS1 before

cache.file <- paste0(cache.dir, "/figA-", Wstyle, ".rds")

dd <- NULL
dd <- if (file.exists(cache.file)) readRDS(cache.file)
if (is.null(dd)) stop("Unable to find the data file", cache.file)

################## Figure
pdf(paste0("fig1.pdf"), width=14, height=5)
	layout(cbind(c(1,1), c(2,2), c(3,5), c(4,6)), widths=c(3,2,2,2))
	par(cex=1, mar=c(4, 4, 3, 1))

	rrr <- do.call(rbind, lapply(dd, function(ddd) sapply(names(phen), function(ppp) whattoconsider(ddd[[ppp]]))))
	prp <- prcomp(rrr, scale.=TRUE)
	plot.PCarrows(prp)
	subpanel("A", line=0)
	
	plot.PCbarplot(prp, mgp=c(2,1,0))
	subpanel("B", line=0)
	
	par(mar=c(3.5,4,1.5,1))
	
	x.mu <- seq(-0.5, 0.5, length.out=11)
	x.sd <- seq( 0, 2, length.out=11)
	x.dd <- seq(0.2, 1, length.out=11)
	x.nn <- seq(2, 12, by=2)
	
	plot(x.mu, sapply(x.mu, function(mm) 100*percentPC(reg.mean=mm)), ylim=c(0,100), xlab=expression(mu[0]), ylab="% variance PC1", mgp=c(2,1,0))
	subpanel("C", line=0.2)
	
	plot(x.sd, sapply(x.sd, function(ss) 100*percentPC(reg.sd=ss)),  ylim=c(0,100), xlab=expression(sigma[0]), ylab="% variance PC1", mgp=c(2,1,0))
	subpanel("D", line=0.2)

	plot(x.dd, sapply(x.dd, function(dd) 100*percentPC(density=dd)), ylim=c(0,100), xlab="Network density", ylab="% variance PC1", mgp=c(2,1,0))
	subpanel("E", line=0.2)
	
	plot(x.nn, sapply(x.nn, function(nn) 100*percentPC(net.size=nn)), ylim=c(0,100), xlab="Network size", ylab="% variance PC1", mgp=c(2,1,0))
	subpanel("F", line=0.2)
	
dev.off()
