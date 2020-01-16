sim.path <- "../src/simul.R"
cache.dir <- "../cache"

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

sim.run <- function(args=default.args, sim.name=NA, force.run=FALSE, nice=TRUE) {
	myargs <- default.args
	myargs[names(args)] <- args # Mistakes in arg names etc. will make the simulation crash later
	if (is.na(sim.name)) 
		sim.name <- paste0("sim", paste(sample(c(letters, LETTERS), 10, replace=TRUE), collapse=""))
	outfile <- paste0(cache.dir, "/", sim.name ,".txt")
	args$outfile <- outfile
	if (force.run || !file.exists(outfile)) {
		comm <- sim.path
		if (nice) comm <- paste0("nice ", sim.path) 
		system(paste(comm, paste0("-", names(args), " ", as.character(unlist(args)), sep="", collapse=" "), sep=" "))
	}
	ans <- read.table(outfile, header=TRUE)
	ans
}
