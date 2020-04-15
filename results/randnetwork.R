

randW.default <- function(net.size, reg.mean, reg.sd, density = 1) {
	# Returns a random matrix of specified size, mean, sd, and density. 
	# Organised in such a way as at least one element in each line is non zero
	numzeros <- min(floor((1-density)*net.size^2), net.size*(net.size-1))
	nonzeros <- (0:(net.size-1))*net.size + sample(1:net.size, net.size, replace=TRUE)
	possiblezeros <- (1:(net.size)^2)[-nonzeros]
	W <- matrix(rnorm(net.size^2, mean=reg.mean, sd=reg.sd), ncol=net.size)
	W[sample(possiblezeros, numzeros, replace=FALSE)] <- 0
	t(W)
}

Wdist.fromWlist <- function(Wlist, epsilon.zero = 0.1) {
	uu <- unlist(Wlist)
	list(
		reg.mean = mean(uu[abs(uu) >= epsilon.zero),
		reg.sd   = sd(uu[abs(uu) >= epsilon.zero),
		density  = mean(abs(uu) < epsilon.zero)
	)
}

Wdist.fromfiles <- function(files, gen = NA, epsilon.zero = 0.1, cols = NA, rows = NA) {
	Wdift.fromWlist(lapply(files, function(ff) {
		ss <- readRDS(ff)
		if (is.na(gen) || !as.character(gen) %in% names(ss)) gen <- rownames(ss)[nrow(ss)]
		W <- ss[[as.character(gen)]]$W
		if (is.na(cols)) cols <- 1:ncol(W)
		if (is.na(rows)) rows <- 1:nrow(W)
		stopifnot(all(cols %in% 1:ncol(W)), all(rows %in% 1:nrow(W)))
		W[rows, cols]
	}), 
	epsilon.zero=epsilon.zero)
} 
