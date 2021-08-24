#!/usr/bin/env Rscript

# Evolutionary trajectories, direct selection

source("./commonsim.R")
source("./terminology.R")
source("./defaults.R")

#################### Options
n.genes          <- default.n
sel.genes        <- default.nsel
s                <- c(rep(default.s, sel.genes), rep(0, n.genes-sel.genes))
W0               <- NA # Random W0

reps             <- default.sim.reps
test.rep         <- default.rob.reps
grad.effect      <- 0.01
N                <- default.N
G                <- 10000
every            <- max(1, round(G/100))
force.run        <- !use.cache

mc.cores         <- default.mc.cores

phen <- c(list(fitness="Fitness"), phen.expression)

col.sim  <- c(oo="black", ie=COL.ENVCAN, le=COL.HOMEO, im=COL.GENCAN, lm=COL.SOM, st=COL.STAB)
col.phen <- c(fitness="black", initenv=COL.ENVCAN, lateenv=COL.HOMEO, initmut=COL.GENCAN, latemut=COL.SOM, stability=COL.STAB)
lty.sim <- c(p=0, m=0, o=0)
pch.sim <- c(p=2, m=6, o=1)
max.points <- 20

######################### Functions
allplots <- function(list.sim, what="fitness", xlim=NULL, ylim=NULL, xlab="Generations", ylab=what, lwd=3, FUN=mean, focal="oo", ...) {
	# Helper function to plot robustness time series
	
	if(is.null(xlim)) 
		xlim <- c(0, as.numeric(rev(names(list.sim[[1]]$mean))[1]))
	allv <- do.call(rbind, lapply(list.sim, function(x) do.call(rbind, lapply(x$mean, function(xx) FUN(xx[[what]])))))
	if(is.null(ylim)) 
		ylim <- c(min(allv), max(allv))
	
	# Reordering the data so that the "focal" traits are plotted at the end and will overlap the others
	if (length(focal) > 0) list.sim <- list.sim[order(sapply(strsplit(names(list.sim), split="\\."), function(x) x[1]) %in% focal)]
	
	plot(NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab="", main=as.expression(phen[what]), col.main=col.phen[what], ...)
	for (nss in names(list.sim)) {
		nnss <- strsplit(nss, split="\\.")[[1]]
		xpp <- seq_along(as.numeric(names(list.sim[[nss]]$mean)))
		xpp <- round(seq(xpp[1], xpp[length(xpp)], length.out=min(max.points, length(xpp))))
		if(any(is.na(sapply(list.sim[[nss]]$mean, function(x) FUN(x[[what]]))[xpp]))) browser()
		lines(as.numeric(names(list.sim[[nss]]$mean))[xpp], sapply(list.sim[[nss]]$mean, function(x) FUN(x[[what]]))[xpp], col=col.sim[nnss[1]], lty=lty.sim[nnss[2]], pch=pch.sim[nnss[2]], type="o", lwd=if(nnss[1] %in% focal) 3 else 1)
	}
}

dynplot <- function(list.sim, what="fitness", focal="oo", xlim=NULL, ylim=NULL, xlab="Generations", ylab=what, lwd=3, FUN=mean, control="oo", ...) {
	if(is.null(xlim)) 
		xlim <- c(0, as.numeric(rev(names(list.sim[[1]]$mean))[1]))
	allv <- do.call(rbind, lapply(list.sim, function(x) do.call(rbind, lapply(x$mean, function(xx) FUN(xx[[what]])))))
	if(is.null(ylim)) 
		ylim <- c(min(allv), max(allv))
	
	plot(NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
	for (nss in names(list.sim)) {
		nnss <- strsplit(nss, split="\\.")[[1]]
		if (! nnss[1] %in% c(focal, control)) next
		xpp <- seq_along(as.numeric(names(list.sim[[nss]]$mean)))
		xpp <- round(seq(xpp[1], xpp[length(xpp)], length.out=min(max.points, length(xpp))))
		lines(
			x   = as.numeric(names(list.sim[[nss]]$mean))[xpp], 
			y   = sapply(list.sim[[nss]]$mean, function(x) FUN(x[[what]]))[xpp], 
			col = col.sim[nnss[1]], 
			lty = lty.sim[nnss[2]], 
			pch = pch.sim[nnss[2]], 
			lwd = if(default.shortcode[what] == focal && nnss[1] == focal) 3 else 1,
			type= "o")
	}
}


################################ Calc
torun <- list(
	oo.o = function() sim.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,0,0,0,0)), 
		reps=reps, series.name="figG-null", force.run=force.run), 
	ie.m = function() sim.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(-grad.effect,0,0,0,0)), 
		reps=reps, series.name="figG-ie-m", force.run=force.run),
	le.m = function() sim.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,-grad.effect,0,0,0)), 
		reps=reps, series.name="figG-le-m", force.run=force.run),
	im.m = function() sim.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,0,-grad.effect,0,0)), 
		reps=reps, series.name="figG-im-m", force.run=force.run),
	lm.m = function() sim.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,0,0,-grad.effect,0)), 
		reps=reps, series.name="figG-lm-m", force.run=force.run),
	st.m = function() sim.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,0,0,0,-grad.effect)), 
		reps=reps, series.name="figG-st-m", force.run=force.run),
	ie.p = function()sim.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(grad.effect,0,0,0,0)), 
		reps=reps, series.name="figG-ie-p", force.run=force.run),
	le.p = function()sim.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,grad.effect,0,0,0)), 
		reps=reps, series.name="figG-le-p", force.run=force.run),
	im.p = function()sim.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,0,grad.effect,0,0)), 
		reps=reps, series.name="figG-im-p", force.run=force.run),
	lm.p = function()sim.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,0,0,grad.effect,0)), 
		reps=reps, series.name="figG-lm-p", force.run=force.run),
	st.p = function() sim.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,0,0,0,grad.effect)), 
		reps=reps, series.name="figG-st-p", force.run=force.run)
)

list.sim <- mclapply(torun, function(ff) ff(), mc.cores=min(length(torun), ceiling(mc.cores/reps)))

########################## Figure

pdf("fig3.pdf", width=12, height=10)
	layout(matrix(1:25, ncol=5))
	par(cex=1, mar=c(0.5, 0.5, 0.2, 0.5), oma=c(5,5,4,3))
	for (selsim in default.shortcode) {
		for (what in names(default.shortcode)) {
			first.col <- (selsim == default.shortcode[1])
			last.col  <- (selsim == default.shortcode[length(default.shortcode)])
			first.row <- (what == names(default.shortcode)[1])
			last.row  <- (what == names(default.shortcode)[length(default.shortcode)])
			
			dynplot(list.sim, what=what, focal=selsim, xlab=if(last.row) "Generation" else "", ylab="", xaxt=if(last.row) "s" else "n", yaxt="n", xpd=NA)
			if (first.col) {
				mtext(default.labels[what], 2, line=1, cex=1.1, las=2)
			}
			if (last.col) {
				axis(4, xpd=NA, las=2)
			}
			if (first.row) {
				selsim.name <- names(default.shortcode)[default.shortcode == selsim]
#~ 				mtext(as.expression(phen.expression[selsim.name]), 3, line=1, col=col.phen[selsim.name], font=2)
				title(as.expression(phen.expression[selsim.name]), line=1, col.main=col.phen[selsim.name], xpd=NA, cex.main=1.1)
			}
			if (first.row && first.col) {
				mtext("Selected robustness component", 3, outer=TRUE, line=2.5, font=2, cex=1.5)
				mtext("Observed robustness component", 2, outer=TRUE, line=3.5, font=2, cex=1.5)
			}
		}
	}
dev.off()
