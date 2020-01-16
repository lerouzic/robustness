sim.path <- "../src/simul.R"
cache.dir <- "../cache"

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
	test.rep=100,
	test.initenv.sd=1,
	test.lateenv.sd=0.1,
	test.initmut.sd=0.1,
	test.lateenv.sd=0.1,
	test.indiv=FALSE)

sim.run <- function(args=default.args, sim.name=NA, force.run=FALSE, nice=TRUE, verbose=FALSE) {
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

acrossrepFUN <- function(simres, FUN=mean, ...) {
	xa <- abind(simres, along=3)
	apply(xa, 1:2, FUN, ...)
}

plotts <- function(mm, sd=NULL, sd.factor=1, colname, xlab="Generation", ylab=colname, col="black", ylim=NULL, xlim=NULL, ...) {
	stopifnot(colname %in% colnames(mm))
	if(dev.cur() == 1) { # no graphical device
		if (is.null(xlim)) xlim <- range(as.numeric(rownames(mm)))
		if (is.null(ylim)) ylim <- range(mm[,colname]) + if(is.null(sd)) 0 else sd.factor*c(-1,1)*max(sd[,colname])
		plot(NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
	}
	lines(x=as.numeric(rownames(mm)), y=mm[,colname], col=col, ...)
	if (!is.null(sd)) {
		arrows(x0=as.numeric(rownames(mm)), y0=mm[,colname]-sd.factor*sd[,colname], y1=mm[,colname]+sd.factor*sd[,colname], lty=1, col=col, length=0)
	}
}
