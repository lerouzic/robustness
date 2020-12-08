
# get the model.M2 function (the genotype-phenotype routine)
curr_dir <- dirname(parent.frame(2)$ofile)   # This kind of things in R always rely on dirty hacks
if (length(curr_dir)==0) curr_dir <- "."
suppressPackageStartupMessages(source(paste(curr_dir, "netw.R", sep="/")))
suppressPackageStartupMessages(source(paste(curr_dir, "robindex.R", sep="/")))

sumlist.u   <- function(ll) Reduce('+', ll)
sumsqlist.u <- function(ll) Reduce('+', lapply(ll, '^', 2)) 
meanlist.u  <- function(ll) sumlist.u(ll)/length(ll)
varlist.u   <- function(ll) sumsqlist.u(ll)/length(ll) - meanlist.u(ll)^2
meanlist    <- function(ll) setNames(lapply(names(ll[[1]]), function(nn) meanlist.u(lapply(ll, "[[", nn))), names(ll[[1]]))
varlist     <- function(ll) setNames(lapply(names(ll[[1]]), function(nn) varlist.u(lapply(ll, "[[", nn))), names(ll[[1]]))

simsel <- function(W0, theta, a=0.2, s=10, grad.rob=rep(0, 5), N=1000, rep=100, G=100, summary.every=1, 
		mut.rate=0.1, som.mut.rate = 0, mut.sd=0.1, initmut.sd=mut.sd, latemut.sd=mut.sd, sim.initenv.sd=0, sim.lateenv.sd=0, 
		initenv.sd=0.1, lateenv.sd=0.1, dev.steps=16, measure=4, plasticity=rep(FALSE, length(theta)), log.robustness=TRUE, mut.correlated=TRUE, ...) {
			
	fitness <- function(x, theta) exp(-sum(s*(x$P-theta)^2) + sum(grad.rob*c(mean(x$initenv), mean(x$lateenv), mean(x$initmut), mean(x$latemut), mean(x$stability))))
	mutate <- function(W) {
		which.mut <- sample(size=1,  which(W != 0)) # Bug if only one W != 0
		W[which.mut] <- rnorm(1, mean=if(mut.correlated) W[which.mut] else 0, sd=mut.sd)
		W
	}
	fmod <- function(x, m) x-m*floor(x/m) 
	renorm <- function(x) {m <- fmod(x,1); x[x<0] <- 1-m[x<0]; x[x>1] <- 1-m[x>1]; x}
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
		plastic <- ifelse(plasticity, runif(1, -0.5, 0.5), 0)
		pop <- lapply(pop, function(i) {
			P <- model.M2(W=i$W, a=a, S0=renorm(rep(a, nrow(i$W)) + plastic + rnorm(nrow(i$W), sd=sim.initenv.sd)), steps=dev.steps, measure=measure)$mean
			if (runif(1) < som.mut.rate) {
				Wmut <- mutate(i$W)
				P <- model.M2(W=Wmut, a=a, S0=P, steps=1, measure=1)$mean
			}
			if (sim.lateenv.sd > 0) {
				P <- model.M2(W=i$W, a=a, S0=renorm(P + rnorm(nrow(i$W), sd=sim.lateenv.sd)), steps=1, measure=1)$mean
			}
			c(i, list(P=P))})

		pop <- lapply(pop, function(i) c(i, list(
			initenv = if (grad.rob[1] != 0 || gg %% summary.every == 0) {
				do.call(robindex.initenv, c(list(W=i$W, a=a, log=log.robustness, rep=rep, 
				env.sd=initenv.sd, dev.steps=dev.steps, measure=measure), dots[names(args(robindex.initenv)) %in% names(dots)]))
			} else { rep(0, nrow(W0)) },
			lateenv = if (grad.rob[2] != 0 || gg %% summary.every == 0) {
				do.call(robindex.lateenv, c(list(W=i$W, a=a, log=log.robustness, rep=rep, 
				env.sd=lateenv.sd, dev.steps=dev.steps, measure=measure), dots[names(args(robindex.lateenv)) %in% names(dots)]))
			} else { rep(0, nrow(W0)) },
			initmut = if (grad.rob[3] != 0 || gg %% summary.every == 0) {
				do.call(robindex.initmut, c(list(W=i$W, a=a, log=log.robustness, rep=rep, 
				mut.sd=initmut.sd, dev.steps=dev.steps, measure=measure), dots[names(args(robindex.initmut)) %in% names(dots)]))
			} else { rep(0, nrow(W0))} ,
			latemut = if (grad.rob[4] != 0 || gg %% summary.every == 0) {
				do.call(robindex.latemut, c(list(W=i$W, a=a, log=log.robustness, rep=rep, 
				mut.sd=latemut.sd, dev.steps=dev.steps, measure=measure), dots[names(args(robindex.latemut)) %in% names(dots)]))
			} else { rep(0, nrow(W0)) },
			stability = if (grad.rob[5] != 0 || gg %% summary.every == 0) {
				do.call(robindex.stability, c(list(W=i$W, a=a, log=log.robustness, dev.steps=dev.steps, measure=measure), 
				dots[names(args(robindex.stability)) %in% names(dots)]))
			} else { rep(0, nrow(W0)) }
			)))
		pop <- lapply(pop, function(i) c(i, list(fitness=fitness(i, theta=renorm(theta+plastic)))))
		if (gg==1 || gg %% summary.every == 0) mpop[[as.character(gg)]] <- meanpop(pop)
		pop <- lapply(pop[sample(seq_along(pop), N, replace=TRUE, prob=sapply(pop, "[", "fitness"))], function(i) list(W=i$W))
	}
	mpop
}
