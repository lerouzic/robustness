#!/usr/bin/env Rscript

# Evolutionary trajectories, direct selection

source("../src/tools.R")
source("./commonsim.R")
source("./terminology.R")
source("./defaults.R")

#################### Options

param <- default

param$s          <- c(rep(param$s, param$nsel), rep(0, param$n-param$nsel))
param$G          <- 10000
param$summary.every<- max(1, round(param$G/100))

W0               <- NA # Random W0
grad.effect      <- 0.01

cache.tag <- "figG"

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

dynplot <- function(list.sim, what="fitness", focal="oo", xlim=NULL, ylim=NULL, xlab="Generations", ylab=what, lwd=3, FUN=mean, control="oo", bg=NA, std.dev=TRUE, relative=FALSE, ...) {
	if(is.null(xlim)) 
		xlim <- c(0, as.numeric(rev(names(list.sim[[1]]$mean))[1]))
	allv  <- do.call(rbind, lapply(list.sim, function(x) do.call(rbind, lapply(x$mean, function(xx) FUN(xx[[what]])))))
	allvv <- do.call(rbind, lapply(list.sim, function(x) do.call(rbind, lapply(x$var,  function(xx) FUN(xx[[what]])))))
	mean.control <- FUN(list.sim[["oo.o"]]$mean[[1]][[what]])
	var.control  <- FUN(list.sim[["oo.o"]]$var [[1]][[what]])
	if (relative) {
		allv  <- allv  - mean.control
		allvv <- allvv - var.control
		allvv[allvv < 0] <- 0 # Negative variances are problematic. 
	}
	if(is.null(ylim)) 
		ylim <- if (std.dev) c(min(allv - sqrt(allvv)), max(allv + sqrt(allvv)))
		        else         c(min(allv), max(allv))
	
	plot(NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
	if (!is.na(bg)) # Not very clean: change the background color of the plotting area
		rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = bg)
		
	jtt <- 0
	for (nss in names(list.sim)) {
		nnss <- strsplit(nss, split="\\.")[[1]]
		if (! nnss[1] %in% c(focal, control)) next

		xpp <- seq_along(as.numeric(names(list.sim[[nss]]$mean)))
		xpp <- round(seq(xpp[1], xpp[length(xpp)], length.out=min(max.points, length(xpp))))
		mean.line <- sapply(list.sim[[nss]]$mean, function(x) FUN(x[[what]]))[xpp]
		if (relative)
			mean.line <- mean.line - mean.control
		lines(
			x   = as.numeric(names(list.sim[[nss]]$mean))[xpp] + jtt, 
			y   = mean.line, 
			col = if (nnss[1] == control) "black" else col.sim[default.shortcode[what]], 
			lty = lty.sim[nnss[2]], 
			pch = pch.sim[nnss[2]], 
			lwd = 1,
			type= "o")
		if (std.dev && nnss[1] != control) {
			var.line <- 
				if (relative) {
					# Here things get more complicated. We want the variance among replicates, relative to the mean at the first generation
					apply(do.call(rbind, lapply(list.sim[[nss]]$full, function(rep) sapply(rep[names(list.sim[[nss]]$mean)], function(gen) FUN(gen[[what]])) - FUN(rep[[1]][[what]]))), 2, var)
				} else {
					sapply(list.sim[[nss]]$var, function(x) FUN(x[[what]]))[xpp]
				}
			arrows(
				x0 = as.numeric(names(list.sim[[nss]]$mean))[xpp] + jtt, 
				y0 = mean.line - sqrt(var.line),
				y1 = mean.line + sqrt(var.line),
				col = makeTransparent(col.sim[default.shortcode[what]], alpha=70), 
				length=0)
			jtt <- jtt + 100
		}
	}
}


################################ Running simulations

torun <- list(
	oo.o = function() sim.run.reps(W0, list(s=param$s, G=param$G, N=param$N, rep=param$rob.reps, summary.every=param$summary.every, grad.rob=c(0,0,0,0,0)), 
		reps=param$sim.reps, series.name=paste0(cache.tag, "-null"), force.run=!param$use.cache), 
	ie.m = function() sim.run.reps(W0, list(s=param$s, G=param$G, N=param$N, rep=param$rob.reps, summary.every=param$summary.every, grad.rob=c(-grad.effect,0,0,0,0)), 
		reps=param$sim.reps, series.name=paste0(cache.tag, "-ie-m"), force.run=!param$use.cache), 
	le.m = function() sim.run.reps(W0, list(s=param$s, G=param$G, N=param$N, rep=param$rob.reps, summary.every=param$summary.every, grad.rob=c(0,-grad.effect,0,0,0)), 
		reps=param$sim.reps, series.name=paste0(cache.tag, "-le-m"), force.run=!param$use.cache), 
	im.m = function() sim.run.reps(W0, list(s=param$s, G=param$G, N=param$N, rep=param$rob.reps, summary.every=param$summary.every, grad.rob=c(0,0,-grad.effect,0,0)), 
		reps=param$sim.reps, series.name=paste0(cache.tag, "-im-m"), force.run=!param$use.cache), 
	lm.m = function() sim.run.reps(W0, list(s=param$s, G=param$G, N=param$N, rep=param$rob.reps, summary.every=param$summary.every, grad.rob=c(0,0,0,-grad.effect,0)), 
		reps=param$sim.reps, series.name=paste0(cache.tag, "-lm-m"), force.run=!param$use.cache), 
	st.m = function() sim.run.reps(W0, list(s=param$s, G=param$G, N=param$N, rep=param$rob.reps, summary.every=param$summary.every, grad.rob=c(0,0,0,0,-grad.effect)), 
		reps=param$sim.reps, series.name=paste0(cache.tag, "-st-m"), force.run=!param$use.cache), 
	ie.p = function()sim.run.reps(W0, list(s=param$s, G=param$G, N=param$N, rep=param$rob.reps, summary.every=param$summary.every, grad.rob=c(grad.effect,0,0,0,0)), 
		reps=param$sim.reps, series.name=paste0(cache.tag, "-ie-p"), force.run=!param$use.cache), 
	le.p = function()sim.run.reps(W0, list(s=param$s, G=param$G, N=param$N, rep=param$rob.reps, summary.every=param$summary.every, grad.rob=c(0,grad.effect,0,0,0)), 
		reps=param$sim.reps, series.name=paste0(cache.tag, "-le-p"), force.run=!param$use.cache), 
	im.p = function()sim.run.reps(W0, list(s=param$s, G=param$G, N=param$N, rep=param$rob.reps, summary.every=param$summary.every, grad.rob=c(0,0,grad.effect,0,0)), 
		reps=param$sim.reps, series.name=paste0(cache.tag, "-im-p"), force.run=!param$use.cache), 
	lm.p = function()sim.run.reps(W0, list(s=param$s, G=param$G, N=param$N, rep=param$rob.reps, summary.every=param$summary.every, grad.rob=c(0,0,0,grad.effect,0)), 
		reps=param$sim.reps, series.name=paste0(cache.tag, "-lm-p"), force.run=!param$use.cache), 
	st.p = function() sim.run.reps(W0, list(s=param$s, G=param$G, N=param$N, rep=param$rob.reps, summary.every=param$summary.every, grad.rob=c(0,0,0,0,grad.effect)), 
		reps=param$sim.reps, series.name=paste0(cache.tag, "-st-p"), force.run=!param$use.cache)
)

list.sim <- mclapply(torun, function(ff) ff(), mc.cores=min(length(torun), ceiling(param$mc.cores/param$sim.reps)))

########################## Figure

pdf("fig3.pdf", width=param$maxfigwidth/param$figscale, height=12/param$figscale, pointsize=param$pointsize)
	layout(matrix(1:25, ncol=5))
	par(cex=1, mar=c(0.5, 0.5, 0.2, 0.5), oma=c(5,5,4,3))
	for (selsim in default.shortcode) {
		for (what in names(default.shortcode)) {
			first.col <- (selsim == default.shortcode[1])
			last.col  <- (selsim == default.shortcode[length(default.shortcode)])
			first.row <- (what == names(default.shortcode)[1])
			last.row  <- (what == names(default.shortcode)[length(default.shortcode)])
			
			dynplot(list.sim, what=what, focal=selsim, xlab=if(last.row) "Generation" else "", ylab="", xaxt=if(last.row) "s" else "n", yaxt="n", xpd=NA, bg=if(selsim != default.shortcode[what]) "gray92" else NA)
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
