
# get the model.M2 function (the genotype-phenotype routine)
library(funr)
curr_dir <- funr::get_script_path()
curr_dir <- paste(curr_dir, "../src", sep="/")
suppressPackageStartupMessages(source(paste(curr_dir, "netw.R", sep="/")))


######## Helper functions
rtruncnorm <- function(N, mean = 0, sd = 1, a = -Inf, b = Inf) {
  if (a > b) stop('Error: Truncation range is empty');
  U <- runif(N, pnorm(a, mean, sd), pnorm(b, mean, sd));
  qnorm(U, mean, sd); }
  

robindex.initenv <- function(W, a, dev.steps, env.sd, rep=1000, FUN=var) {
	ans <- replicate(rep,  
		model.M2(W, a, S0=rtruncnorm(nrow(W), mean=a, sd=env.sd, 0, 1) , steps=dev.steps, measure=1)$mean)
	apply(ans, 1, FUN)
}

robindex.lateenv <- function(W, a, dev.steps, env.sd, rep=1000, FUN=var) {
	ref <- model.M2(W, a, S0=rep(a, ncol(W)) , steps=dev.steps, measure=1)$mean
	ans <- replicate(rep, 
		model.M2(W, a, S0=rtruncnorm(nrow(W), mean=ref, sd=env.sd, 0, 1) , steps=1, measure=1)$mean)
	apply(ans, 1, FUN)
}

robindex.initmut <- function(W, a, dev.steps, mut.sd, mut.correlated=FALSE, nbmut=1, rep=1000, FUN=var) {
	ans <- replicate(rep,  
		{
			myW <- W
			nbmut <- min(nbmut, sum(W != 0))
			which.mut <- sample(size=nbmut,  which(W != 0), replace=FALSE) # Bug if only one W != 0
			mm <- if(mut.correlated) W[which.mut] else 0
			myW[which.mut] <- rnorm(nbmut, mean=mm, sd=mut.sd)
			model.M2(myW, a, S0=rep(a, ncol(W)), steps=dev.steps, measure=1)$mean
		})
	apply(ans, 1, FUN)
}

robindex.latemut <- function(W, a, dev.steps, mut.sd, mut.correlated=FALSE, nbmut=1, rep=1000, FUN=var) {
	ref <- model.M2(W, a, S0=rep(a, ncol(W)) , steps=dev.steps, measure=1)$mean
	ans <- replicate(rep,  
		{
			myW <- W
			nbmut <- min(nbmut, sum(W != 0))
			which.mut <- sample(size=nbmut,  which(W != 0), replace=FALSE) # Bug if only one W != 0
			mm <- if(mut.correlated) W[which.mut] else 0
			myW[which.mut] <- rnorm(nbmut, mean=mm, sd=mut.sd)
			model.M2(myW, a, S0=ref, steps=1, measure=1)$mean
		})
	apply(ans, 1, FUN)
}

robindex.stability <- function(W, a, dev.steps) {
	ref <- model.M2(W, a, S0=rep(a, ncol(W)) , steps=dev.steps, measure=1)$mean
	onemore <- model.M2(W, a=a, S0=ref, steps=1, measure=1)$mean
	(ref-onemore)^2
}
