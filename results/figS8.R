source("./defaults.R")
source("./randnetwork.R")
source("../src/netw.R")

a             <- default.a
dev.steps     <- default.dev.steps

mut.sd        <- 0.5 #default.sim.mutsd
mut.correlated<- default.mut.correlated

rand.density  <- 0.5 #default.density
rand.mean     <- default.rand.mean
rand.sd       <- default.rand.sd
net.size      <- default.n


# Duplicated from fig J
mutate <- function(W, mut.sd) {
	which.mut <- sample(size=1,  which(W != 0)) # Bug if only one W != 0
	W[which.mut] <- rnorm(1, mean=if(mut.correlated) W[which.mut] else 0, sd=mut.sd)
	W
}

n.mut <- c(1,2,5)
replicates <- 200

set.seed(0123456) # Making sure to get a reproducible example

W <- randW(net.size, rand.mean, rand.sd, rand.density)

pdf("figS8.pdf", height=3, width=3*length(n.mut))
layout(t(seq_along(n.mut)))
for (nn in n.mut) {
	plot(NULL, xlab="Expression gene 1", ylab="Expression gene 2", xlim=0:1, ylim=0:1, main=paste0(nn, " mutation", if(nn>1) "s" else ""))
	for (r in 1:replicates) {
		myW <- W
		for (mm in 1:nn) myW <- mutate(myW, mut.sd)
		mod <- model.M2 (W=myW, a=a, steps=dev.steps)$mean
		points(mod[1], mod[2], pch=1, col="gray")
	}
	mod.ref <- model.M2 (W=W, a=a, steps=dev.steps)$mean
	points(mod.ref[1], mod.ref[2], pch=1, col="blue", cex=3)
}
dev.off()
