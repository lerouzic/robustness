#!/usr/bin/env Rscript

library(parallel)
mc.cores <- min(64, detectCores()-1)

source("./terminology.R")
source("./defaults.R")
source("./studycases.R")
source("../src/robindex.R")

use.cache <- TRUE

phen <- phen.expression

grid.size      <- 101

a              <- default.a
dev.steps      <- default.dev.steps
measure        <- default.dev.measure
log            <- default.log.robustness

rob.reps       <- 10000
rob.initenv.sd <- default.initenv.sd
rob.lateenv.sd <- default.lateenv.sd
rob.initmut.sd <- default.initmut.sd
rob.latemut.sd <- default.latemut.sd

cache.dir <- "../cache"
cache.file <- paste(cache.dir, "figC.rds", sep="/")
if (!dir.exists(cache.dir)) dir.create(cache.dir)

difftarget <- function(res, target) {
    df <- t(res[,grep("mean.", colnames(res))])-target
    apply(abs(df), 2, sum)
}

plotres <- function(res, crit="mean", stud=NULL, mask=NULL, contour=FALSE, mx = 0.2) {
    z <- res[,grep(colnames(res), pattern=crit)[1]] # take only the first gene
    z[z>mx] <- mx
    if (!is.null(mask)) z[mask] <- NA
    
    if (crit %in% names(phen)) main <- phen[crit] else main <- crit
    
    image(x=ww1, y=ww2, z=matrix(z, nrow=sqrt(nrow(res))), main=main, xlab=expression(W[11]), ylab=expression(W[21]), col=heat.colors(128)[1:100])
    if (contour) contour(x=ww1, y=ww2, z=matrix(z, nrow=sqrt(nrow(res))), add=TRUE) 
    if (!is.null(stud))
        invisible(sapply(rownames(stud), function(rn) text(x=stud[rn,1], y=stud[rn,2], rn, col="blue")))
}

ww1 <- seq(-1.5, 2, length.out=grid.size)
ww2 <- seq(-1, 4, length.out=grid.size)
# Pattern 
#   a   NA
#   b   NA
# (should not be important)


if (!use.cache || !file.exists(cache.file)) {
	res <- expand.grid(ww1, ww2)
	
	res <- cbind(res, t(sapply(1:nrow(res), function(i) {
	        w <- targetW(W=cbind(unlist(res[i,]), rep(NA, network.size)), target=target, a=a)
	        return( c(w)[(1+network.size*(network.size-1)):(network.size*network.size)])
	    })))
	colnames(res) <- paste("W", outer(1:network.size, 1:network.size, paste, sep="."), sep=".")
	
	res <- cbind(res, do.call(rbind, mclapply(1:nrow(res), function(i) {
	        w <- matrix(unlist(res[i,]), nrow=network.size)
	        c(mean=model.M2(w, a, steps=dev.steps)$mean,
	          initenv=robindex.initenv(w, a, dev.steps, measure, rob.initenv.sd, rep=rob.reps, log=log),
	          lateenv=robindex.lateenv(w, a, dev.steps, measure, rob.initenv.sd, rep=rob.reps, log=log),
	          initmut=robindex.initmut(w, a, dev.steps, measure, rob.initmut.sd, rep=rob.reps, log=log),
	          latemut=robindex.latemut(w, a, dev.steps, measure, rob.latemut.sd, rep=rob.reps, log=log),
	          stability=robindex.stability(w, a, dev.steps, log=log))
	    }, mc.cores=mc.cores)))
	saveRDS(res, file=cache.file)
} else {
	res <- readRDS(cache.file)
}

pdf("figC.pdf", width=12, height=8)
	layout(rbind(1:3, c(4:5, 0)))
	mm <- difftarget(res, target) > 0.15
	for (ppp in names(phen))
		plotres(res, ppp, stud, mask=mm, contour=TRUE)
dev.off()

