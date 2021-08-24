
# get the model.M2 function (the genotype-phenotype routine)
library(funr)
curr_dir <- funr::get_script_path()
curr_dir <- paste(curr_dir, "../src", sep="/")
suppressPackageStartupMessages(source(paste(curr_dir, "netw.R", sep="/")))


min.log <- -40 # Always embarrassing to deal with NaNs....
mylog <- function(x) ifelse(x < exp(min.log),min.log, log(x))

######## Helper functions
rtruncnorm <- function(N, mean = 0, sd = 1, a = -Inf, b = Inf) {
  if (a > b) stop('Error: Truncation range is empty');
  U <- runif(N, pnorm(a, mean, sd), pnorm(b, mean, sd));
  qnorm(U, mean, sd); }
  

conditional <- function(M, constr) {
	# Conditional evolvability, modified from Hansen & Houle 2008
	if(!is.numeric(constr)) constr <- match(constr, rownames(M))
	stopifnot(sum(is.na(constr)) == 0, all(constr %in% seq_len(ncol(M))))
	nconstr <- seq_len(ncol(M))[!seq_len(ncol(M)) %in% constr]
	My  <- M[nconstr,nconstr, drop=FALSE]
	Myx <- M[nconstr, constr, drop=FALSE]
	Mxy <- M[constr, nconstr, drop=FALSE]
	Mx  <- M[constr, constr, drop=FALSE]
	My - Myx %*% solve(Mx) %*% Mxy
}

robindex.initenv <- function(W, a, dev.steps, measure=min(4, round(dev.steps/5)), env.sd, rep=1000, FUN=var, log=FALSE) {
	ans <- replicate(rep,  
		model.M2(W, a, S0=rtruncnorm(nrow(W), mean=a, sd=env.sd, 0, 1) , steps=dev.steps, measure=measure)$mean)
	transf <- if(log) function(x) mylog(x) else identity   # wierd syntax, just to solve the non-elegant problem of calling a local variable "log"
	transf(apply(ans, 1, FUN))
}

robindex.lateenv <- function(W, a, dev.steps, measure=min(4, round(dev.steps/5)), env.sd, rep=1000, FUN=var, log=FALSE) {
	ref <- model.M2(W, a, S0=rep(a, ncol(W)) , steps=dev.steps, measure=measure)$mean
	ans <- replicate(rep, 
		model.M2(W, a, S0=rtruncnorm(nrow(W), mean=ref, sd=env.sd, 0, 1) , steps=1, measure=1)$mean)
	transf <- if(log) function(x) mylog(x) else identity
	transf(apply(ans, 1, FUN))
}

robindex.initmut <- function(W, a, dev.steps, measure=min(4, round(dev.steps/5)), mut.sd, mut.correlated=TRUE, nbmut=1, rep=1000, FUN=var, log=FALSE) {
	ans <- replicate(rep,  
		{
			myW <- W
			nbmut <- min(nbmut, sum(W != 0))
			which.mut <- sample(size=nbmut,  which(W != 0), replace=FALSE) # Bug if only one W != 0
			mm <- if(mut.correlated) W[which.mut] else 0
			myW[which.mut] <- rnorm(nbmut, mean=mm, sd=mut.sd)
			model.M2(myW, a, S0=rep(a, ncol(W)), steps=dev.steps, measure=measure)$mean
		})
	transf <- if(log) function(x) mylog(x) else identity
	transf(apply(ans, 1, FUN))
}

robindex.latemut <- function(W, a, dev.steps, measure=min(4, round(dev.steps/5)), mut.sd, mut.correlated=TRUE, nbmut=1, rep=1000, FUN=var, log=FALSE) {
	ref <- model.M2(W, a, S0=rep(a, ncol(W)) , steps=dev.steps, measure=measure)$mean
	ans <- replicate(rep,  
		{
			myW <- W
			nbmut <- min(nbmut, sum(W != 0))
			which.mut <- sample(size=nbmut,  which(W != 0), replace=FALSE) # Bug if only one W != 0
			mm <- if(mut.correlated) W[which.mut] else 0
			myW[which.mut] <- rnorm(nbmut, mean=mm, sd=mut.sd)
			model.M2(myW, a, S0=ref, steps=1, measure=1)$mean
		})
	transf <- if(log) function(x) mylog(x) else identity
	transf(apply(ans, 1, FUN))
}

robindex.stability <- function(W, a, dev.steps, measure=min(4, round(dev.steps/5)), log=FALSE) {
	ref <- model.M2(W, a, S0=rep(a, ncol(W)) , steps=dev.steps, measure=measure)$mean
	onemore <- model.M2(W, a=a, S0=ref, steps=1, measure=1)$mean
	transf <- if(log) function(x) mylog(x) else identity
	transf((ref-onemore)^2)
}

robindex.Mmatrix <- function(W, a, dev.steps, mut.sd=0.1, mut.correlated=FALSE, test.initmut.sd=mut.sd, test.latemut.sd=mut.sd, nbmut=1, test.initenv.sd=1, test.lateenv.sd=0.1, test.rep=100, rep=1000, log.robustness=FALSE, include.expr=FALSE) {
	all <- replicate(rep, {
		myW <- W
		nbmut <- min(nbmut, sum(W != 0))
		which.mut <- sample(size=nbmut,  which(W != 0), replace=FALSE)
		mm <- if(mut.correlated) W[which.mut] else 0
		myW[which.mut] <- rnorm(nbmut, mean=mm, sd=mut.sd)
		c(
		expr=if(include.expr) model.M2(W=myW, a=a, S0=rep(a, ncol(W)) , steps=dev.steps, measure=min(4, round(dev.steps/5)))$mean,
		initenv=mean(robindex.initenv(W=myW, a=a, dev.steps=dev.steps, env.sd=test.initenv.sd, rep=test.rep, log=log.robustness)),
		lateenv=mean(robindex.lateenv(W=myW, a=a, dev.steps=dev.steps, env.sd=test.lateenv.sd, rep=test.rep, log=log.robustness)),
		initmut=mean(robindex.initmut(W=myW, a=a, dev.steps=dev.steps, mut.sd=test.initmut.sd, mut.correlated=mut.correlated, rep=test.rep, log=log.robustness)),
		latemut=mean(robindex.latemut(W=myW, a=a, dev.steps=dev.steps, mut.sd=test.latemut.sd, mut.correlated=mut.correlated, rep=test.rep, log=log.robustness)),
		stability=mean(robindex.stability(W=myW, a=a, dev.steps=dev.steps, log=log.robustness)))
	})
	list(mean=rowMeans(all), vcov=var(t(all), na.rm=TRUE)/nbmut)
}

robindex.Mmatrix.outfile <- function(out, gen=NA, ...) {
	# out is the data.frame corresponding to an output file
	# Warning: it probably does not make sense to call this function on genotypes averaged across simulations!
	stopifnot(is.data.frame(out))
	tags <- c("initenv", "lateenv", "initmut", "latemut", "stability")
	if (is.na(gen) || !(gen %in% as.numeric(rownames(out)))) gen <- as.numeric(rownames(out))
	ans <- lapply(gen, function(gg) {
		mm <- sapply(tags, function(tt) mean(unlist(out[as.character(gg), grepl(colnames(out), pattern=paste("robustness", tt, sep="."))])))
		ww <- unlist(out[as.character(gg), grepl(colnames(out), pattern="genotype.mean")])
		ww <- matrix(ww, ncol=sqrt(length(ww)))
		mmat <- robindex.Mmatrix(W=ww, ...)
		colnames(mmat$vcov) <- rownames(mmat$vcov) <- tags
		list(mean=mm, vcov=mmat$vcov)
	})
	names(ans) <- as.character(gen)
	ans
}

robindex.Gmatrix.outfile <- function(out, gen=NA, ...) {
	# out is the data.frame corresponding to an output file
	stopifnot(is.data.frame(out))
	tags <- c("initenv", "lateenv", "initmut", "latemut", "stability")
	if (is.na(gen) || !(gen %in% as.numeric(rownames(out)))) gen <- as.numeric(rownames(out))
	ans <- lapply(gen, function(gg) {
		mm <- sapply(tags, function(tt) mean(unlist(out[as.character(gg), grepl(colnames(out), pattern=paste("robustness", tt, sep="."))])))
		gmat <- unlist(out[as.character(gg), grepl(colnames(out), pattern="robustness.G")])
		stopifnot(length(mm) == sqrt(length(gmat)))
		gmat <- matrix(gmat, ncol=length(mm))
		colnames(gmat) <- rownames(gmat) <- tags
		list(mean=mm, vcov=gmat)
	})
	names(ans) <- as.character(gen)
	ans
}
