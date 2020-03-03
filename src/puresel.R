

meanlist <- function(ll) lapply(names(ll[[1]]), function(nn) 

puresel <- function(W0, theta, a=0.2, s=10, grad.rob=rep(0, 5), N=1000, rep=100, G=100, mut.rate=0.1, mut.sd=0.1, log=TRUE, ...) {
	fitness <- function(x) exp(-s*sum(x$P-theta)^2 + grad.rob*c(x$initenv, x$lateenv, x$initmut, x$latemut, x$stability))
	meanpop <- function(pop) meanlist(pop)
	
	stopifnot(is.matrix(W0), nrow(W0) > 0, nrow(W0) == ncol(W0))
	stopifnot(length(theta) == ncol(W0))
	stopifnot(length(s) == 1 || length(s) == length(theta))
	stopifnot(length(grad.rob) == 5)
	
	dots <- as.list(...)
	
	pop <- lapply(1:N, list(W=W0))
	
	for (gg in 1:G) {
		pop <- lapply(pop, function(i) c(i, list(P=model.M2(W=i$W, a=a, steps=20)$mean)))
		pop <- lapply(pop, function(i) c(i, list(
			initenv=do.call(robindex.initenv, c(list(W=i$W, a=a, log=log), dots[names(args(robindex.initenv)) %in% names(dots)])),
			lateenv=do.call(robindex.lateenv, c(list(W=i$W, a=a, log=log), dots[names(args(robindex.lateenv)) %in% names(dots)])),
			initmut=do.call(robindex.initmut, c(list(W=i$W, a=a, log=log), dots[names(args(robindex.initmut)) %in% names(dots)])),
			latemut=do.call(robindex.latemut, c(list(W=i$W, a=, log=loga), dots[names(args(robindex.latemut)) %in% names(dots)])),
			stability=do.call(robindex.stability, c(list(W=i$W, a=a, log=log), dots[names(args(robindex.stability)) %in% names(dots)])))))
		pop <- lapply(pop, function(i) c(i, list(fitness=fitness(i)))
		
	}
	
	
}
