#!/usr/bin/env Rscript

# Ensures that the last PCs are real (not due to stochastic noise)

source("./terminology.R")
source("./defaults.R")

source("../src/randnetwork.R")
source("../src/robindex.R")
source("../src/pca.R")
source("../src/tools.R")

######################## Options

param <- default

cache.tag <- "figH"

param$rob.reps  <- 1000
net.reps        <- 10000     # Number of simulated networks, should match figure S2

allreps         <- round(10^(seq(2, 4, length.out=7)))
allrobs         <- round(10^(seq(1, 5, length.out=9)))

cols <- 1:5

########################### Running simulations

# Technically, this sounds pretty useless: it would be way more efficient to run the max number of replicates and sample them for lower
# counts. Yet, the current code is simpler, and avoids correlations among samples. 

cache.file.reps <- paste0(param$cache.dir, "/", cache.tag, "-reps.rds")
cache.file.robs <- paste0(param$cache.dir, "/", cache.tag, "-robs.rds")

res.reps <- if (param$use.cache && file.exists(cache.file.reps)) {
				readRDS(cache.file.reps)
			} else {
				ans <- lapply(allreps, function(reps) {
					Wlist <- replicate(reps, randW(net.size=param$n, reg.mean=param$rand.mean, reg.sd=param$rand.sd, density=param$density), simplify=FALSE)
					eigenV(Wlist, rob.reps=param$rob.reps, param=param, summary.FUN=param$summary.FUN)
				})
				saveRDS(ans, file=cache.file.reps, version=2)
				ans
			}

res.robs <- if (param$use.cache && file.exists(cache.file.robs)) {
				readRDS(cache.file.robs)
			} else {
				ans <- lapply(allrobs, function(robs) {
					Wlist <- replicate(net.reps, randW(net.size=param$n, reg.mean=param$rand.mean, reg.sd=param$rand.sd, density=param$density), simplify=FALSE)
					eigenV(Wlist, rob.reps=robs, param=param, summary.FUN=param$summary.FUN)
				})
				saveRDS(ans, file=cache.file.robs, version=2)
				ans
			}

############################ Figure

pdf(paste0("figS3.pdf"), width=12/param$figscale, height=6/param$figscale, pointsize=param$pointsize)
	layout(t(1:2))
	
	plot(NULL, xlim=c(0.4,1)*range(allreps), ylim=c(1e-3,1), log="xy", xlab="Number of simulated networks", ylab="Proportion variance explained", yaxt="n", xaxt="n")
	for (i in 1:5)
		lines(allreps, sapply(res.reps, function(r) r[i]), col=cols[i], type="o", pch=16)
	text(allreps[1], res.reps[[1]], col=cols, pos=2, paste0("PC", seq_along(cols)))
	axis(2, at=c(0.001,0.01,0.1,1), labels=c("0.1%", "1%", "10%", "100%"))
	axis(1, at=c(100, 1000, 10000))
	subpanel("A", line=1)
	
	plot(NULL, xlim=c(0.2,1)*range(allrobs), ylim=c(1e-3,1), log="xy", xlab="Number of robustness tests", ylab="Proportion variance explained", yaxt="n", xaxt="n")
	for (i in 1:5)
		lines(allrobs, sapply(res.robs, function(r) r[i]), col=cols[i], type="o", pch=16)
	text(allrobs[1], res.robs[[1]], col=cols, pos=2, paste0("PC", seq_along(cols)))	
	axis(2, at=c(0.001,0.01,0.1,1), labels=c("0.1%", "1%", "10%", "100%"))
	axis(1, at=c(10, 1000, 100000), labels=c("10", "1000", "100000"))
	subpanel("B", line=1)
	
dev.off()

