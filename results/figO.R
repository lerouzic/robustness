source("./randnetwork.R")
source("./defaults.R")
source("../src/netw.R")

library(parallel)

cache.dir <- "../cache"
use.cache <- TRUE

dev.steps   <- default.dev.steps
dev.measure <- default.dev.measure
a           <- default.a
net.size    <- default.n

density     <- default.density
reps        <- 1000

mc.cores    <- default.mc.cores

rand.mean   <- default.rand.mean
rand.sd     <- default.rand.sd

mm.range <- c(min(-0.5, 0.9*rand.mean), max(0.5, 1.1*rand.mean))
ss.range <- c(0, max(0.5, 1.1*rand.sd))

cache.file <- paste0(cache.dir, "/figO.rds")
if (!dir.exists(cache.dir)) dir.create(cache.dir)	

resolution <- 41

mm.all <- seq(mm.range[1], mm.range[2], length.out=resolution)
ss.all <- seq(ss.range[1], ss.range[2], length.out=resolution)

ii <- NULL
ii <- if (use.cache && file.exists(cache.file)) readRDS(cache.file)

if (is.null(ii)) {
	ii <- outer(mm.all, ss.all, function(mm, ss) mapply(mm, ss, FUN=function(mmm, sss) {
			ee <- unlist(mclapply(1:reps, function(i)  model.M2(randW(net.size, mmm, sss, density), a, steps=dev.steps, measure=dev.measure)$mean, mc.cores=mc.cores))
			ks.test(ee, "punif")$statistic
		}))
	saveRDS(ii, file=cache.file)
}


pdf("figO.pdf", width=10, height=5)
	layout(t(1:2))
	
	image(x=mm.all, y=ss.all, z=ii, xlab="Mean W", ylab="Std.dev. W", col = hcl.colors(1024, "YlOrRd", rev = TRUE) ,main="Kolmogorov-Smirnoff test statistic")
	contour(x=mm.all, y=ss.all, z=ii, add=TRUE)
	points(rand.mean, rand.sd, pch=1, col="red", cex=2, lwd=3)
	
	ee.ref <- unlist(mclapply(1:(100*reps), function(i)  model.M2(randW(net.size, rand.mean, rand.sd, density), a, steps=dev.steps, measure=dev.measure)$mean, mc.cores=mc.cores))
	hist(ee.ref, breaks=50, xlab="Gene expression", freq=FALSE, main="Distribution of gene expressions")
dev.off()
