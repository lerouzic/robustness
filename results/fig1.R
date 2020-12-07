#!/usr/bin/env Rscript

# Runs a PCA among robustness indexes and plots the results

source("./terminology.R")
source("./defaults.R")

#################### Options

Wstyle <- "random" # Possible: "random", "evolved" ,"randevol"
whattoconsider <- function(x) mean(x) # the average index for all genes

phen <- phen.expression      # from terminology.R    
phen.pos.ref <- "initmut"    # this guy will always be positive (so that the figure is reproducible)

################### Functions
myplot.prcomp <- function(pr, labels=default.labels, cols=default.cols) {
    par(xpd=NA)
    nPC <- length(pr$sdev)
    layout(t(1:2))
    
    # In prcomp, the sign of the PCs is arbitrary. To ensure reproducibility, PCs will be oriented in such a was that phen.pos.ref > 0
    pr$rotation <- t(t(pr$rotation)*sign(pr$rotation[phen.pos.ref,]))
    
    plot(NULL, xlab="", ylab="", ylim=c(0.75,nPC+0.25), xlim=range(pr$rotation), yaxt="n", bty="n")
    arrows(x0=min(pr$rotation), y0=1:nPC, x1=max(pr$rotation), code=3, length=0.2)
    for (i in 1:nPC) text(x=pr$rotation[,i], y=nPC-i+1+0.1*(0:(nPC-1)), labels[rownames(pr$rotation)], col=cols[rownames(pr$rotation)], pos=3, cex=0.8, font=2)
    lines(x=rep(0,2), y=c(1,nPC), lty=2)
    axis(2, at=1:nPC, labels=paste0("PC", nPC:1))
    
    bb <- barplot(100*rev(pr$sdev^2/(sum(pr$sdev^2))), horiz=TRUE, ylim=c(1, 2*nPC), ylab="", xlab="% variance", width=1, space=1)
    arrows(x0=seq(20, 80, 20), y0=1, y1=max(bb), col="gray", lty=3, length=0)
}

#################### Calc
source("./figS1.R") # This script uses dataset figA. Run figS1 before

cache.file <- paste0(cache.dir, "/figA-", Wstyle, ".rds")

dd <- NULL
dd <- if (file.exists(cache.file)) readRDS(cache.file)
if (is.null(dd)) stop("Unable to find the data file", cache.file)

################## Figure
pdf(paste0("fig1.pdf"), width=7, height=5)
	rrr <- do.call(rbind, lapply(dd, function(ddd) sapply(names(phen), function(ppp) whattoconsider(ddd[[ppp]]))))
	prp <- prcomp(rrr, scale.=TRUE)
	myplot.prcomp(prp)
dev.off()
