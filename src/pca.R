source("../src/robindex.R")

eigenV <- function(Wlist, rob.reps, param, summary.FUN=mean, index=c("initenv", "lateenv", "initmut", "latemut", "stability")) {
	dd <- mclapply(Wlist, function(W) {
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
			test.reps       = rob.reps, 
			log.robustness  = param$log.robustness)
	}, mc.cores=param$mc.cores) 

	rrr <- do.call(rbind, lapply(dd, function(ddd) sapply(index, function(ppp) summary.FUN(ddd[[ppp]]))))
	prp <- prcomp(rrr, scale.=TRUE)
	prp$sdev^2/(sum(prp$sdev^2))
}
