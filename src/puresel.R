


# get the model.M2 function (the genotype-phenotype routine)
curr_dir <- dirname(parent.frame(2)$ofile)   # This kind of things in R always rely on dirty hacks
if (length(curr_dir)==0) curr_dir <- "."
suppressPackageStartupMessages(source(paste(curr_dir, "netw.R", sep="/")))
suppressPackageStartupMessages(source(paste(curr_dir, "robindex.R", sep="/")))

sumlist.u   <- function(ll) Reduce('+', ll)
sumsqlist.u <- function(ll) Reduce(function(x,y) x^2+y^2, ll)
meanlist.u  <- function(ll) sumlist.u(ll)/length(ll)
varlist.u   <- function(ll) sumsqlist.u(ll)/length(ll) - meanlist.u(ll)^2
meanlist    <- function(ll) setNames(lapply(names(ll[[1]]), function(nn) meanlist.u(lapply(ll, "[[", nn))), names(ll[[1]]))
varlist     <- function(ll) setNames(lapply(names(ll[[1]]), function(nn) varlist.u(lapply(ll, "[[", nn))), names(ll[[1]]))

puresel <- function(W0, theta, a=0.2, s=10, grad.rob=rep(0, 5), N=1000, rep=100, G=100, summary.every=1, 
		mut.rate=0.1, mut.sd=0.1, initmut.sd=mut.sd, latemut.sd=mut.sd, initenv.sd=0.1, lateenv.sd=0.1, 
		dev.steps=16, log.robustness=TRUE, mut.correlated=TRUE, ...) {
			
	fitness <- function(x) exp(-sum(s*(x$P-theta)^2) + sum(grad.rob*c(mean(x$initenv), mean(x$lateenv), mean(x$initmut), mean(x$latemut), mean(x$stability))))
	mutate <- function(W) {
		which.mut <- sample(size=1,  which(W != 0)) # Bug if only one W != 0
		W[which.mut] <- rnorm(1, mean=if(mut.correlated) W[which.mut] else 0, sd=mut.sd)
		W
	}
	meanpop <- function(pop) meanlist(pop)
	
	stopifnot(is.matrix(W0), nrow(W0) > 0, nrow(W0) == ncol(W0))
	stopifnot(length(theta) == 1 || length(theta) == ncol(W0))
	stopifnot(length(s) == 1 || length(s) == ncol(W0))
	stopifnot(length(grad.rob) == 5)
	
	if (length(theta)==1) theta <- rep(theta, ncol(W0))
	if (length(s) == 1) s <- rep(s, ncol(W0))
	
	dots <- list(...)
	
	pop <- lapply(1:N, function(i) list(W=W0))
	mpop <- list()
	
	for (gg in 1:G) {
		pop <- lapply(pop, function(i) { if (runif(1) < mut.rate) i$W <- mutate(i$W); i })
		pop <- lapply(pop, function(i) c(i, list(P=model.M2(W=i$W, a=a, steps=20)$mean)))
		pop <- lapply(pop, function(i) c(i, list(
			initenv=do.call(robindex.initenv, c(list(W=i$W, a=a, log=log.robustness, rep=rep, env.sd=initenv.sd, dev.steps=dev.steps), dots[names(args(robindex.initenv)) %in% names(dots)])),
			lateenv=do.call(robindex.lateenv, c(list(W=i$W, a=a, log=log.robustness, rep=rep, env.sd=lateenv.sd, dev.steps=dev.steps), dots[names(args(robindex.lateenv)) %in% names(dots)])),
			initmut=do.call(robindex.initmut, c(list(W=i$W, a=a, log=log.robustness, rep=rep, mut.sd=initmut.sd, dev.steps=dev.steps), dots[names(args(robindex.initmut)) %in% names(dots)])),
			latemut=do.call(robindex.latemut, c(list(W=i$W, a=a, log=log.robustness, rep=rep, mut.sd=latemut.sd, dev.steps=dev.steps), dots[names(args(robindex.latemut)) %in% names(dots)])),
			stability=do.call(robindex.stability, c(list(W=i$W, a=a, log=log.robustness, dev.steps=dev.steps), dots[names(args(robindex.stability)) %in% names(dots)])))))
		pop <- lapply(pop, function(i) c(i, list(fitness=fitness(i))))
		if (gg %% summary.every == 0) mpop[[as.character(gg)]] <- meanpop(pop)
		pop <- lapply(pop[sample(seq_along(pop), N, replace=TRUE, prob=sapply(pop, "[", "fitness"))], function(i) list(W=i$W))
	}
	mpop
}
