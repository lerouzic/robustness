#!/usr/bin/env Rscript

# Runs a PCA among robustness indexes and plots the results

source("./terminology.R")
source("./defaults.R")

source("../src/robindex.R") 
source("../src/randnetwork.R")
source("../src/pca.R")
source("../src/tools.R")

source("./figS2.R") # This script uses dataset figA. Run figS2 before


#################### Options

param <- default

param$rob.reps     <- 1000
net.reps           <- 10000

cache.tag <- "figA"

cols.PC <- 1:5
show.PC <- 1:3
  
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

percentPC <- function(reps=net.reps, net.size=default$n, density=default$density, reg.mean=default$rand.mean, reg.sd=default$rand.sd, PC=1:net.size) {
	cache.str  <- paste (cache.tag, net.size, density, reg.mean, reg.sd, sep="-")
	cache.file <- paste0(param$cache.dir, "/", cache.str, ".rds")
	all.Wrob <-  
		if (file.exists(cache.file)) readRDS(cache.file)
		else {
			allW <- lapply(1:reps, function(i) randW(net.size=net.size, reg.mean=reg.mean, reg.sd=reg.sd, density=density))
			ans  <- eigenV(allW, param$rob.reps, param, summary.FUN=param$summary.FUN)
			saveRDS(ans, file=cache.file, version=2)
			ans
		}
	all.Wrob
}

#################### Importing the dataset from figS2

cache.file <- paste0(param$cache.dir, "/", cache.tag, "-main.rds")

dd.fS2 <- if (file.exists(cache.file)) {
			readRDS(cache.file)
		  } else {
			stop("Unable to find the data file", cache.file)
		  }
		  
################### Figure

pdf(paste0("fig1.pdf"), width=14, height=5)

	layout(cbind(c(1,1), c(2,2), c(3,5), c(4,6)), widths=c(3,2,2,2))
	par(cex=1, mar=c(4, 4, 3, 1))

	rrr <- do.call(rbind, lapply(dd.fS2, function(ddd) sapply(names(phen.expression), function(ppp) param$summary.FUN(ddd[[ppp]]))))
	prp <- prcomp(rrr, scale.=TRUE)
	plot.PCarrows(prp)
	subpanel("A", line=0)
	
	plot.PCbarplot(prp, mgp=c(2,1,0))
	subpanel("B", line=0)
	
	par(mar=c(3.5,4,1.5,1))
	
	x.mu <- seq(-0.5, 0.5, length.out=11)
	x.sd <- seq( 0,   2,   length.out=11)
	x.dd <- seq( 0.2, 1,   length.out=11)
	x.nn <- seq( 2,   12,   by=2)
	
	plot(NULL, xlim=range(x.mu), ylim=c(0,100), xlab=expression(mu[0]), ylab="% variance PC", mgp=c(2,1,0))
	for (ppc in show.PC)
		lines(x.mu, sapply(x.mu, function(mm) 100*percentPC(reg.mean=mm)[ppc]), col=cols.PC[ppc])
	legend("left", lty=1, col=cols.PC[show.PC], legend=paste0("PC", show.PC), bty="n")
	subpanel("C", line=0.2)
	
	plot(NULL, xlim=range(x.sd), ylim=c(0,100), xlab=expression(sigma[0]), ylab="% variance PC", mgp=c(2,1,0))
	for (ppc in show.PC)
		lines(x.sd, sapply(x.sd, function(ss) 100*percentPC(reg.sd=ss)[ppc]), col=cols.PC[ppc])
	subpanel("D", line=0.2)

	plot(NULL, xlim=range(x.dd), ylim=c(0,100), xlab="Network density", ylab="% variance PC", mgp=c(2,1,0))
	for (ppc in show.PC)
		lines(x.dd, sapply(x.dd, function(dd) 100*percentPC(density=dd)[ppc]), col=cols.PC[ppc])	
	subpanel("E", line=0.2)
	
	plot(NULL, xlim=range(x.nn), ylim=c(0,100), xlab="Network size", ylab="% variance PC", mgp=c(2,1,0))
	for (ppc in show.PC)
		lines(x.nn, sapply(x.nn, function(nn) 100*percentPC(net.size=nn)[ppc]), col=cols.PC[ppc])	
	subpanel("F", line=0.2)
	
dev.off()
