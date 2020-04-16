

randW <- function(net.size, reg.mean, reg.sd, density = 1) {
	# Returns a random matrix of specified size, mean, sd, and density. 
	# Organised in such a way as at least one element in each line is non zero
	numzeros <- min(floor((1-density)*net.size^2), net.size*(net.size-1))
	nonzeros <- (0:(net.size-1))*net.size + sample(1:net.size, net.size, replace=TRUE)
	possiblezeros <- (1:(net.size)^2)[-nonzeros]
	W <- matrix(rnorm(net.size^2, mean=reg.mean, sd=reg.sd), ncol=net.size)
	W[sample(possiblezeros, numzeros, replace=FALSE)] <- 0
	t(W)
}

Wdist.fromWlist <- function(Wlist, epsilon.zero = 0.01) {
	uu <- unlist(Wlist)
	list(
		mean = mean(uu[abs(uu) >= epsilon.zero]),
		sd   = sd(uu[abs(uu) >= epsilon.zero]),
		density  = mean(abs(uu) > epsilon.zero)
	)
}

Wlist.fromfiles <- function(files, gen = NA, cols = NA, rows = NA) {
	lapply(files, function(ff) {
		ss <- readRDS(ff)
		if (is.na(gen) || !as.character(gen) %in% names(ss)) gen <- names(ss)[length(ss)]
		W <- ss[[as.character(gen)]]$W
		if (is.na(cols)) cols <- 1:ncol(W)
		if (is.na(rows)) rows <- 1:nrow(W)
		stopifnot(all(cols %in% 1:ncol(W)), all(rows %in% 1:nrow(W)))
		W[rows, cols]
	})
} 

Wdist.fromfiles <- function(files, gen = NA, epsilon.zero = 0.1, cols = NA, rows = NA) {
	Wdist.fromWlist(Wlist.fromfiles(files=files, gen=gen, cols=cols, rows=rows), epsilon.zero=epsilon.zero)
} 
