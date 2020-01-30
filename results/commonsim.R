sim.path <- "../src/simul.R"
cache.dir <- "../cache"

library(parallel)
library(abind)

default.args <- list(
	a=0.2,
	r=0.5,
	n=6, 
	init.sd=0.0001,
	mut.correlated=TRUE,
	dev.steps=20,
	G=100,
	summary.every=10,
	mc.cores=1,
	N=100,
	mut.rate=0.01,
	mut.sd=0.1,
	som.rate=0,
	som.sd=0,
	initenv.sd=0,
	lateenv.sd=0,
	theta=0.5,
	s=10,
	ss=0,
	test.rep=100,
	test.initenv.sd=1,
	test.lateenv.sd=0.1,
	test.initmut.sd=0.1,
	test.latemut.sd=0.1,
	test.indiv=TRUE,
	log.robustness=TRUE)
	
sim.run.single <- function(args=default.args, sim.name=NA, force.run=FALSE, nice=TRUE, verbose=FALSE) {
	myargs <- default.args
	myargs[names(args)] <- args # Mistakes in arg names etc. will make the simulation crash later
	if (is.na(sim.name)) 
		sim.name <- paste0("sim", paste(sample(c(letters, LETTERS), 10, replace=TRUE), collapse=""))
	outfile <- paste0(cache.dir, "/", sim.name ,".txt")
	myargs$outfile <- outfile
	if (force.run || !file.exists(outfile)) {
		comm <- sim.path
		if (nice) comm <- paste0("nice ", sim.path) 
		comm <-paste(comm, paste0("-", names(myargs), " ", as.character(unlist(myargs)), sep="", collapse=" "), sep=" ") 
		if (verbose) cat(comm, "\n")
		system(comm)
	}
	ans <- read.table(outfile, header=TRUE)
	ans
}

sim.run.reps <- function(args=NULL, reps=10, series.name="sim", force.run=FALSE, nice=TRUE, verbose=FALSE, mc.cores=detectCores()-1) {
	ans <- mclapply(1:reps, function(r) {
		ss <- sim.run.single(args, paste0(series.name, "-", r), force.run=force.run, nice=nice, verbose=verbose)
	}, mc.cores=mc.cores)

	ans.mean <- acrossrepFUN(ans, mean)
	ans.sd   <- acrossrepFUN(ans, sd)
	list(full=ans, mean=as.data.frame(ans.mean), sd=as.data.frame(ans.sd))
}

acrossrepFUN <- function(simres, FUN=mean, ...) {
	xa <- abind(simres, along=3)
	apply(xa, 1:2, FUN, ...)
}

lighten.color <- function(color, factor=0.3){
    col <- col2rgb(color)
    col <- col+(255-col)*factor
    col <- rgb(t(col), maxColorValue=255)
    col
}


plotts <- function(mm, sd=NULL, sd.factor=1, colname, xlab="Generation", ylab=colname, col="black", ylim=NULL, xlim=NULL, add=FALSE, type="l", ...) {
	stopifnot(colname %in% colnames(mm))
	if(dev.cur() == 1 || !add) { # no graphical device
		if (is.null(xlim)) xlim <- range(as.numeric(rownames(mm)))
		if (is.null(ylim)) ylim <- range(mm[,colname]) + if(is.null(sd)) 0 else sd.factor*c(-1,1)*max(sd[,colname])
		plot(NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
	}
	if (!is.null(sd)) {
		arrows(x0=as.numeric(rownames(mm)), y0=mm[,colname]-sd.factor*sd[,colname], y1=mm[,colname]+sd.factor*sd[,colname], lty=1, col=lighten.color(col), length=0)
	}
	lines(x=as.numeric(rownames(mm)), y=mm[,colname], col=col, type=type, ...)
}

plotmat <- function(listmat, arg1=1, arg2=2, col="black", col.end=col, sd.factor=2, var.thresh=1e-25, add=FALSE, ...) {
	library(ellipse)
	if (!add) {
		xlim <- c(min(sapply(listmat, function(x) x$mean[arg1] - 1.2*sd.factor*sqrt(x$vcov[arg1, arg1]))), max(sapply(listmat, function(x) x$mean[arg1] + 1.2*sd.factor*sqrt(x$vcov[arg1, arg1])), na.rm=TRUE))
		ylim <- c(min(sapply(listmat, function(x) x$mean[arg2] - 1.2*sd.factor*sqrt(x$vcov[arg2, arg2]))), max(sapply(listmat, function(x) x$mean[arg2] + 1.2*sd.factor*sqrt(x$vcov[arg2, arg2])), na.rm=TRUE))
		plot(NULL, xlim=xlim, ylim=ylim, xlab=arg1, ylab=arg2, ...)
	}
	for (i in seq_along(listmat)) {
		if (listmat[[i]]$vcov[arg1,arg1] < var.thresh) { 
			listmat[[i]]$vcov[arg1,arg1] <- var.thresh
			listmat[[i]]$vcov[arg1,arg2] <- listmat[[i]]$vcov[arg2,arg1] <- 0
		}
		if (listmat[[i]]$vcov[arg2,arg2] < var.thresh) { 
			listmat[[i]]$vcov[arg2,arg2] <- var.thresh
			listmat[[i]]$vcov[arg1,arg2] <- listmat[[i]]$vcov[arg2,arg1] <- 0
		}
		lines(ellipse(cov2cor(listmat[[i]]$vcov[c(arg1, arg2),][,c(arg1,arg2)]), 
			scale=sqrt(diag(listmat[[i]]$vcov)[c(arg1, arg2)]), level=2*pnorm(sd.factor)-1,
			centre=listmat[[i]]$mean[c(arg1, arg2)]), col=rgb(colorRamp(c(col, col.end))(i/length(listmat))/255))
	}
}
