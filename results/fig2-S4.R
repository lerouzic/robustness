#!/usr/bin/env Rscript

# For convenience, two figures (fig 2 and fig S4) are plotted at the same time. 

source("./terminology.R")
source("./defaults.R")
source("./studycases.R") # defines network.size and target

source("../src/robindex.R")
source("../src/tools.R")

################## Options

param <- default

param$rob.reps <- 10000

cache.tag         <- "figC"
grid.size         <- 101
difftarget.thresh <- 0.15


################# Functions

difftarget <- function(res, target) {
	# Difference between the network output and the target (simply the sum of the absolute values of the differences)
    df <- t(res[,grep("mean.", colnames(res))])-target
    apply(abs(df), 2, sum)
}

whyitfails <- function(W, a, dev.steps, measure, target) {
	# 0: OK
	# 1: Still evolving
	# 2: Stuck to the border
	# 3: Limit cycle
	# -1: Unknown
	evol.thresh <- 0.002
	border.thresh <- 0.02
	mm <- model.M2(W=W, a=a, steps=dev.steps, measure=measure, full=TRUE)
	Smes <- mm$full[,(dev.steps-measure+2):(dev.steps+1)]
	if (sum(abs(mm$mean-target)) <= difftarget.thresh) 
		return (0)
	if (any(apply(Smes, 1, function(x) length(unique(sign(diff(x))))==1 && all(abs(diff(x)) >= evol.thresh))))
		return(1)
	if (any(apply(Smes, 1, function(x) all(x < border.thresh) || all(x > 1 - border.thresh))))
		return(2)
	if (any(apply(Smes, 1, function(x) { lb <- x < border.thresh; lu <- x > 1 - border.thresh; return((any(lb) && !all(lb)) || (any(lu) && !all(lu)))})))
		return(3)	
	return(-1)
}

plotres <- function(res, crit="mean", stud=NULL, mask=NULL, contour=FALSE, mx = 0.2) {
	# Helper function for the 2D plot
    z <- apply(res[,grep(colnames(res), pattern=crit)], 1,  param$summary.FUN)# take only the first gene
    z[z>mx] <- mx
    if (!is.null(mask)) z[mask] <- NA
    
    if (crit %in% names(phen.expression)) main <- as.expression(phen.expression[crit]) else main <- crit
    
    image(x=ww1, y=ww2, z=matrix(z, nrow=sqrt(nrow(res))), main=main, xlab=expression(W[11]), ylab=expression(W[21]), col=heat.colors(128)[1:100])
    if (contour) contour(x=ww1, y=ww2, z=matrix(z, nrow=sqrt(nrow(res))), add=TRUE, col="gray45") 
    if (!is.null(stud))
        invisible(sapply(rownames(stud), function(rn) { 
			points(x=stud[rn,1], y=stud[rn,2], pch=16, col=makeTransparent("white", 95) ,cex=3)
			text(x=stud[rn,1], y=stud[rn,2], rn, col="blue") 
		}))
}

##################### Simulations

cache.file <- paste0(param$cache.dir, "/", cache.tag, ".rds")

ww1 <- seq(-1.5, 2, length.out=grid.size)
ww2 <- seq(-1,   4, length.out=grid.size)
# Pattern 
#   a   NA
#   b   NA
# (should not be important)

if (!param$use.cache || !file.exists(cache.file)) {
	res <- expand.grid(ww1, ww2)
	
	res <- cbind(res, t(sapply(1:nrow(res), function(i) {
	        w <- targetW(W=cbind(unlist(res[i,]), rep(NA, network.size)), target=target, a=param$a)
	        return( c(w)[(1+network.size*(network.size-1)):(network.size*network.size)])
	    })))
	colnames(res) <- paste("W", outer(1:network.size, 1:network.size, paste, sep="."), sep=".")
	
	res <- cbind(res, do.call(rbind, mclapply(1:nrow(res), function(i) {
	        w <- matrix(unlist(res[i,]), nrow=network.size)
	        robindex.Wmatrix(w, param$a, param$dev.steps, param$dev.measure, 
				mut.correlated=param$mut.correlated, 
				test.initmut.sd=param$initmut.sd, test.latemut.sd=param$latemut.sd, test.initenv.sd=param$initenv.sd, test.lateenv.sd=param$lateenv.sd, 
				test.reps=param$rob.reps, log.robustness=param$log.robustness)
	    }, mc.cores=param$mc.cores)))

	res <- cbind(res, WIF=apply(res, 1, function(x) {
		w <-  matrix(unlist(x[1:(network.size^2)]), nrow=network.size)
		whyitfails(w, param$a, param$dev.steps, param$dev.measure, target)
	}))
	saveRDS(res, file=cache.file, version=2)
} else {
	res <- readRDS(cache.file)
}


######################## Figures

pdf("fig2.pdf", width=param$maxfigwidth/param$figscale, height=10/param$figscale, pointsize=param$pointsize)
	layout(rbind(1:3, c(4:5, 0)))
	par(cex=1, mar=c(3,3,2,0.5), mgp=c(2,1,0))
	mm <- difftarget(res, target) > difftarget.thresh
	for (ppp in names(phen.expression))
		plotres(res, ppp, stud, mask=mm, contour=TRUE)
dev.off()

pdf("figS4.pdf", width=6.5/param$figscale, height=6/param$figscale, pointsize=param$pointsize)
	zz <- matrix(res$WIF, nrow=sqrt(nrow(res)))
	zz[zz==-1] <- 2 # When stuck to the border, this is still an alternative equilibrium
	cc <- c('On target' = "white", 'Still changing'="yellow", 'Alternative eq.'="gray", 'Large osc.'="red")
	image(x=ww1, y=ww2, z=zz, xlab=expression(W[11]), ylab=expression(W[21]), col=cc)
	legend("topleft", inset=c(0,-0.1), pch=15, col=cc, legend=names(cc), horiz=TRUE, xpd=TRUE, cex=0.8, pt.cex=2, bty="n")
dev.off()
