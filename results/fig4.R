#!/usr/bin/env Rscript

source("./commonsim.R")
source("./terminology.R")
source("./defaults.R")
source("../src/robindex.R")
cache.dir <- "../cache"

library(parallel)
library(abind)

mc.cores <- default.mc.cores

use.cache <- TRUE

n.genes           <- default.n
sel.genes         <- default.nsel    
s                 <- c(rep(default.s, sel.genes), rep(0, n.genes-sel.genes))
W0                <- NA
reps              <- default.sim.reps
test.rep          <- default.rob.reps
grad.effect       <- 0.01
N                 <- default.N 
n                 <- default.n
a                 <- default.a
G                 <- 10000
every             <- round(G/100)
force.run         <- !use.cache

pch.sim <- c(o=1, m=6, p=2, mm=10, pp=11, mp=9, pm=7)

############################ Helper (plotting) functions
plot.rel <- function(list.sim, what, G=NULL, pch=NULL, xlab="Relative response", ...) {
	# list.sim: the list of simulation results
	# what: full length names of robustness indicators to plot (can be a vector)
	
	if (is.null(G)) 
		G <- rev(names(list.sim[[1]]$mean))[1] # Default: last generation
		
	if (is.null(pch))
		pch <- pch.sim[c("mm","mp","pm","pp")]

	plot(NULL, xlim=c(-2,1), ylim=c(0.5,0.8 + length(what)), xlab=xlab, ylab="", yaxt="n", bty="n", ...)
	axis(2, at=seq_along(what), tick=FALSE, labels=as.expression(phen.expression[what]), las=2, mgp=c(3,0,0))
	abline(v=0, lty=3, col="darkgray")
	
	for (i in seq_along(what)) {
		ref.o <- mean(list.sim[["oo.o"]]$mean[[G]][[what[i]]])
		ref.m <- mean(list.sim[[paste0(default.shortcode[what[i]], ".m")]]$mean[[G]][[what[i]]])
		ref.p <- mean(list.sim[[paste0(default.shortcode[what[i]], ".p")]]$mean[[G]][[what[i]]])

		arrows(
			x0=(ref.m-ref.o)/(ref.p-ref.m), 
			x1=(ref.p-ref.o)/(ref.p-ref.m), 
			y0=i, 
			code=3, angle=90, length=0.05, col=default.cols[what[i]], lwd=3)

		for (j in seq_along(what)) {
			if (j == i) next
			ss.name <- paste0(default.shortcode[what[i]], ".", default.shortcode[what[j]])
			ss <- sapply(names(pch), function(nn) mean(list.sim[[paste0(ss.name, ".", nn)]]$mean[[G]][[what[i]]]))
			points((ss-ref.o)/(ref.p-ref.m), rep(i+0.1*j, length(ss)), pch=pch[names(ss)], col=default.cols[what[j]])
		}
	}
	legend("topright", pch=pch, 
		legend=c("target - & corr -     ", "target - & corr +     ", "target + & corr -", "target + & corr +    "), 
		horiz=TRUE, bty="n", xpd=NA, cex=0.8)
}

plot.ts <- function(list.sim, what, w1, w2, xlab="Generation", ylab=as.expression(phen.expression[what]), ylim=NULL, lwd=2, ...) {
	# list.sim: the list of simulation results
	# what: full length name of the robustness indicator
	# w1  : full length name of the first selected robustness indicator (target trait), should be "what" most of the time
	# w2  : fill length name of the second selected rob indicator (corr trait).
	
	Gmax <- as.numeric(rev(names(list.sim[[1]]$mean))[1])
	xx <- unique(round(seq(from=1, to=length(list.sim[[1]]$mean), length.out=15)))
	
	sw1 <- default.shortcode[w1]
	sw2 <- default.shortcode[w2]
	
	ref.o <- sapply(list.sim[["oo.o"]]$mean, function(x) mean(x[[what]]))
	ref.m <- sapply(list.sim[[paste0(sw1, ".m")]]$mean, function(x) mean(x[[what]]))
	ref.p <- sapply(list.sim[[paste0(sw1, ".p")]]$mean, function(x) mean(x[[what]]))
	
	sim.mm <- sapply(list.sim[[paste0(sw1, ".", sw2, ".mm")]]$mean, function(x) mean(x[[what]]))
	sim.mp <- sapply(list.sim[[paste0(sw1, ".", sw2, ".mp")]]$mean, function(x) mean(x[[what]]))
	sim.pm <- sapply(list.sim[[paste0(sw1, ".", sw2, ".pm")]]$mean, function(x) mean(x[[what]]))
	sim.pp <- sapply(list.sim[[paste0(sw1, ".", sw2, ".pp")]]$mean, function(x) mean(x[[what]]))
	
	if (is.null(ylim))
		ylim <- range(c(ref.m, ref.p, sim.mm, sim.pp))
	
	plot(NULL, xlim=c(1, Gmax), ylim=ylim, xlab=xlab, ylab=ylab, ...)
	lines(as.numeric(names(ref.o))[xx], ref.o[xx], pch=pch.sim["o"], lwd=lwd)
	lines(as.numeric(names(ref.p))[xx], ref.p[xx], lty=1, col=default.cols[w1], lwd=lwd) #pch=pch.sim["p"])
	lines(as.numeric(names(ref.m))[xx], ref.m[xx], lty=1, col=default.cols[w1], lwd=lwd) #pch=pch.sim["m"])
	
	points(as.numeric(names(sim.mm))[xx], sim.mm[xx], pch=pch.sim["mm"], col=default.cols[w2])
	points(as.numeric(names(sim.mp))[xx], sim.mp[xx], pch=pch.sim["mp"], col=default.cols[w2])
	points(as.numeric(names(sim.pm))[xx], sim.pm[xx], pch=pch.sim["pm"], col=default.cols[w2])
	points(as.numeric(names(sim.pp))[xx], sim.pp[xx], pch=pch.sim["pp"], col=default.cols[w2])
}

########################### Simulations

# Computes gradient verctors for one and two selected traits
gradvec1 <- function(i, grd) c(rep(0, i-1), grd, rep(0, length(default.shortcode)-i))
gradvec2 <- function(i1, i2, grd1, grd2) c(rep(0, i1-1), grd1, rep(0, i2-i1-1), grd2, rep(0, length(default.shortcode)-i2))

# Make a list of simulation names and gradients
torun <- list(
	oo.o = list(series.name="figG-null", grad.rob=c(0,0,0,0,0)))
for (i in seq_along(default.shortcode)) {
	nm <- default.shortcode[i]
	torun[[paste(nm, "m", sep=".")]] <- list(series.name=paste("figG", nm, "m", sep="-"), grad.rob=gradvec1(i, -grad.effect))
	torun[[paste(nm, "p", sep=".")]] <- list(series.name=paste("figG", nm, "p", sep="-"), grad.rob=gradvec1(i,  grad.effect))
}
for (i1 in 1:(length(default.shortcode)-1))
	for (i2 in (i1+1):length(default.shortcode)) {
		nm1 <- default.shortcode[i1]
		nm2 <- default.shortcode[i2]
		torun[[paste(nm1, nm2, "mm", sep=".")]] <- list(series.name=paste("figI", nm1, nm2, "mm", sep="-"), grad.rob=gradvec2(i1, i2, -grad.effect, -grad.effect))
		torun[[paste(nm1, nm2, "mp", sep=".")]] <- list(series.name=paste("figI", nm1, nm2, "mp", sep="-"), grad.rob=gradvec2(i1, i2, -grad.effect,  grad.effect))
		torun[[paste(nm1, nm2, "pm", sep=".")]] <- list(series.name=paste("figI", nm1, nm2, "pm", sep="-"), grad.rob=gradvec2(i1, i2,  grad.effect, -grad.effect))
		torun[[paste(nm1, nm2, "pp", sep=".")]] <- list(series.name=paste("figI", nm1, nm2, "pp", sep="-"), grad.rob=gradvec2(i1, i2,  grad.effect,  grad.effect))
	}

list.sim <- mclapply(torun, function(ff) 
	sim.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=ff$grad.rob), reps=reps, series.name=ff$series.name, force.run=force.run, mc.cores=min(mc.cores, reps)), 
	mc.cores=ceiling(mc.cores/reps))

# It is more convenient to have all combinations in the list.sim variable (should not take more memory)
for (i in 1:(length(default.shortcode)-1))
	for (j in (i+1):length(default.shortcode)) {
		nc <- paste0(default.shortcode[i], ".", default.shortcode[j], ".")
		inc <- paste0(default.shortcode[j], ".", default.shortcode[i], ".")
		list.sim[[paste0(inc, "pp")]] <- list.sim[[paste0(nc, "pp")]]
		list.sim[[paste0(inc, "mp")]] <- list.sim[[paste0(nc, "pm")]]
		list.sim[[paste0(inc, "pm")]] <- list.sim[[paste0(nc, "mp")]]
		list.sim[[paste0(inc, "mm")]] <- list.sim[[paste0(nc, "mm")]]
	}

####################### Figure

pdf("fig4.pdf", width=8, height=5)
	layout(cbind(c(1,1,1),c(2,3,4)), widths=c(2,1))
	par(mar=c(4, 10, 1, 1), cex=0.75, oma=c(2,0,0,0))
	plot.rel(list.sim, names(default.shortcode), xlab="")
	mtext(side=1, "Relative response", xpd=NA, line=3, cex=0.75)
	
	# The example panels are not encoded programatically, manual adjustments required
	arrows(x0=0.5, y0=5.4, x1=1.2, y1=5.4, lwd=2, xpd=NA, length=0.1, col=default.cols["latemut"])
	arrows(x0=0.5, y0=3.2, x1=1.3, y1=3.25, lwd=2, xpd=NA, length=0.1, col=default.cols["lateenv"])
	arrows(x0=0.7, y0=1.5, x1=1.3, y1=1.3, lwd=2, xpd=NA, length=0.1, col=default.cols["stability"])
	
	par(mar=c(2, 4, 1, 1))
	plot.ts(list.sim, "stability", "stability", "latemut", xlab="", mgp=c(2,1,0))
	plot.ts(list.sim, "initmut", "initmut", "lateenv", xlab="", mgp=c(2,1,0), ylim=c(-11,-6))
	plot.ts(list.sim, "initenv", "initenv", "stability", xlab="", mgp=c(2,1,0), ylim=c(-38,-6))
	mtext(side=1, "Generations", xpd=NA, line=2.5, cex=0.75)
	
dev.off()



