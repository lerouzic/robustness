#!/usr/bin/env Rscript

source("./commonsim.R")
source("./terminology.R")
source("./defaults.R")

source("../src/robindex.R")


param <- default

param$s           <- c(rep(param$s, param$nsel), rep(0, param$n-param$nsel))
param$G           <- 10000
param$summary.every<- round(param$G/100)

W0                <- NA
grad.effect       <- 0.01

M.reps            <- 100
M.nbmut           <- 5

max.points <- 12

cache.tag.uni <- "figG"
cache.tag.bi  <- "figI"
cache.tag.mut <- "figP"

col.sim  <- c(oo="black", ie=COL.ENVCAN, le=COL.HOMEO, im=COL.GENCAN, lm=COL.SOM, st=COL.STAB)
pch.sim <- c(o=1, m=6, p=2, mm=10, pp=11, mp=9, pm=7)

rob.pairs <- list(
	c("ie","le"),
	c("ie","im"),
	c("ie","lm"),
	c("ie","st"),
	c("le","im"),
	c("le","lm"),
	c("le","st"),
	c("im","lm"),
	c("im","st"),
	c("lm","st")
)


#################### Functions

# Computes gradient verctors for one and two selected traits
# (code puplicated from fig 4)
gradvec1 <- function(i, grd) c(rep(0, i-1), grd, rep(0, length(default.shortcode)-i))
gradvec2 <- function(i1, i2, grd1, grd2) c(rep(0, i1-1), grd1, rep(0, i2-i1-1), grd2, rep(0, length(default.shortcode)-i2))

meancor <- function(Mlist, i1, i2, nm=NULL, mc.cores=1) { # param$mc.cores) {
	larr <- mclapply(Mlist, function(Mrep) {
		unlist(sapply(Mrep, function(Mgen) cov2cor(Mgen$vcov)[i1,i2]))
	}, mc.cores = mc.cores)

	if(is.null(nm)) 
		nm <- unique(unlist(lapply(rev(larr), names)))
	larr <- lapply(larr, function(x) setNames(x[nm], nm))
	
	arr <- do.call(cbind, larr)
	rowMeans(arr, na.rm=TRUE)
}

################### Running simulations

# Make a list of simulation names and gradients
torun <- list(
	oo.o = list(series.name=paste0(cache.tag.uni, "-null"), grad.rob=c(0,0,0,0,0)))
for (nm in unique(unlist(rob.pairs))) {
	torun[[paste(nm, "m", sep=".")]] <- list(series.name=paste0(cache.tag.uni, "-", nm, "-m"), grad.rob=gradvec1(which(default.shortcode==nm), -grad.effect))
	torun[[paste(nm, "p", sep=".")]] <- list(series.name=paste0(cache.tag.uni, "-", nm, "-p"), grad.rob=gradvec1(which(default.shortcode==nm),  grad.effect))
}
for (icomp in rob.pairs) {
	nm1 <- icomp[1]
	nm2 <- icomp[2]
	i1  <- which(default.shortcode==nm1)
	i2  <- which(default.shortcode==nm2)
	stopifnot(i1 < i2)
	torun[[paste(nm1, nm2, "mm", sep=".")]] <- list(series.name=paste(cache.tag.bi, nm1, nm2, "mm", sep="-"), grad.rob=gradvec2(i1, i2, -grad.effect, -grad.effect))
	torun[[paste(nm1, nm2, "mp", sep=".")]] <- list(series.name=paste(cache.tag.bi, nm1, nm2, "mp", sep="-"), grad.rob=gradvec2(i1, i2, -grad.effect,  grad.effect))
	torun[[paste(nm1, nm2, "pm", sep=".")]] <- list(series.name=paste(cache.tag.bi, nm1, nm2, "pm", sep="-"), grad.rob=gradvec2(i1, i2,  grad.effect, -grad.effect))
	torun[[paste(nm1, nm2, "pp", sep=".")]] <- list(series.name=paste(cache.tag.bi, nm1, nm2, "pp", sep="-"), grad.rob=gradvec2(i1, i2,  grad.effect,  grad.effect))
}

# Running simulations
#  in theory, if fig4 has been run before, no need for new simulations, read from the cache instead. 
list.sim <- mclapply(torun, function(ff) 
	sim.run.reps(W0, list(
			s            = param$s, 
			G            = param$G, 
			N            = param$N, 
			rep          = param$rob.reps, 
			summary.every= param$summary.every, 
			grad.rob     = ff$grad.rob), 
		reps         = param$sim.reps, 
		series.name  = ff$series.name, 
		force.run    = !param$use.cache, 
		mc.cores     = max(1,floor(param$sim.reps/param$mc.cores))), 
	mc.cores = ceiling(param$mc.cores/param$sim.reps))

# Computing M matrices
list.M <- mclapply(setNames(nm=names(list.sim)), function(name.sim.series) {
			cache.file <- paste0(param$cache.dir, "/", cache.tag.mut, "-", name.sim.series, ".rds")
			if (param$use.cache && file.exists(cache.file)) {
				return(readRDS(cache.file))
			} else {
				ans <- mclapply(list.sim[[name.sim.series]]$full, function(simrep) {
					lapply(simrep, function(gen) {
						robindex.Mmatrix(
							W              = gen$W,
							a              = param$a, 
							dev.steps      = param$dev.steps, 
							mut.sd         = param$sim.mutsd, 
							mut.correlated = param$mut.correlated, 
							test.initmut.sd= param$initmut.sd, 
							test.latemut.sd= param$latemut.sd, 
							nbmut          = M.nbmut, 
							test.initenv.sd= param$initenv.sd, 
							test.lateenv.sd= param$lateenv.sd, 
							test.rep       = param$rob.reps, 
							rep            = M.reps, 
							log.robustness = param$log.robustness, 
							include.expr   = TRUE)
					})
				}, mc.cores=max(1,floor(param$mc.cores/length(list.sim)))
				)
				saveRDS(ans, cache.file, version=2)
				ans
			} 
		}, mc.cores=param$mc.cores)


################### Figure

cairo_pdf("fig5.pdf", width=12/param$figscale, height=10/param$figscale, pointsize=param$pointsize)
	lm <- matrix(0, ncol=4, nrow=4)
	lm[lower.tri(lm, diag=TRUE)] <- 1:10
	lm[1,4] <- 11   # top-right panel
	layout(lm)
	
	par(cex=1, mar=c(1, 1, 0.5, 0.5), oma=c(5, 4, 1, 0.5))
	
	for (icomp in seq_along(rob.pairs)) {
		nm1 <- rob.pairs[[icomp]][1]
		nm2 <- rob.pairs[[icomp]][2]
		i1  <- which(default.shortcode==nm1)
		i2  <- which(default.shortcode==nm2)
		xx <- as.numeric(unique(unlist(lapply(rev(list.M[["oo.o"]]), names))))
		xx <- xx[unique(round(seq(1, length(xx), length.out=max.points)))]
		plot(NULL, xlim=c(0, param$G), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n")
		# , main=substitute(r(r1, r2), list(r1=as.list(default.labels[i1])[[1]], r2=as.list(default.labels[i2])[[1]]))
		if (nm1 == "ie") {
			axis(2, at=seq(0, 1, by=0.25))
			mtext(as.expression(phen.expression[names(i2)]), 2, line=3.5, xpd=NA, col=default.cols[names(i2)])
			if (nm2 == "le")
				mtext("Mutational correlation", side=2, line=1.3, outer=TRUE)
		}
		if (nm2 == "st") {
			axis(1, at=c(0, 5000, 10000))
			mtext(as.expression(phen.expression[names(i1)]), 1, line=3.5, xpd=NA, col=default.cols[names(i1)])
			if (nm1 == "ie")
				mtext("Generation", side=1, line=1.3, outer=TRUE)
		}
		
		points(xx, meancor(list.M[["oo.o"]], names(i1), names(i2), nm=as.character(xx)), col=col.sim["oo"], pch=pch.sim["o"])
		
		points(xx, meancor(list.M[[paste0(nm1, ".p")]], names(i1), names(i2), nm=as.character(xx)), col=col.sim[nm1], pch=pch.sim["p"])
		points(xx, meancor(list.M[[paste0(nm1, ".m")]], names(i1), names(i2), nm=as.character(xx)), col=col.sim[nm1], pch=pch.sim["m"])
		points(xx, meancor(list.M[[paste0(nm2, ".p")]], names(i1), names(i2), nm=as.character(xx)), col=col.sim[nm2], pch=pch.sim["p"])
		points(xx, meancor(list.M[[paste0(nm2, ".m")]], names(i1), names(i2), nm=as.character(xx)), col=col.sim[nm2], pch=pch.sim["m"])
	}
	
	# top-right panel
	plot(NULL, xlim=c(-0.2,0.2), ylim=c(-0.15, 0.3), xaxt="n", yaxt="n", 
		xlab=expression(Delta*r*" selection for lower "*rho), 
		ylab=expression(Delta*r*" selection for larger "*rho), 
		xpd=NA)
	axis(1, at=c(-0.2, 0, 0.2))
	axis(2, at=c(0, 0.2))
	rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "bisque")
	abline(h=0, v=0, col="darkgray")
	
	for (lnm1 in names(default.shortcode)) {
		for (lnm2 in names(default.shortcode)) {
			if (lnm2 == lnm1) next
			ww1 <- which(names(default.shortcode) == lnm1)
			ww2 <- which(names(default.shortcode) == lnm2)
			if (ww2 > ww1) ww2 <- ww2 - 1
			yy.m <- diff(meancor(list.M[[paste0(default.shortcode[lnm1], ".m")]], lnm1, lnm2, nm=as.character(xx[c(1,length(xx))])))
			yy.p <- diff(meancor(list.M[[paste0(default.shortcode[lnm1], ".p")]], lnm1, lnm2, nm=as.character(xx[c(1,length(xx))])))
			
			points(yy.m, yy.p, pch=21, col=col.sim[default.shortcode[lnm2]], bg=col.sim[default.shortcode[lnm1]])
		}
	}
dev.off()


cairo_pdf("figS9.pdf", width=10/param$figscale, height=10/param$figscale, pointsize=param$pointsize)
	lm <- matrix(0, ncol=4, nrow=4)
	lm[lower.tri(lm, diag=TRUE)] <- 1:10
	lm[1,3] <- 11   # caption
	lm[2,4] <- 12   # color scale
	layout(lm)
	
	par(cex=1, mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
	
	for (icomp in seq_along(rob.pairs)) {
		nm1 <- rob.pairs[[icomp]][1]
		nm2 <- rob.pairs[[icomp]][2]
		i1  <- which(default.shortcode==nm1)
		i2  <- which(default.shortcode==nm2)
		xx <- as.numeric(unique(unlist(lapply(rev(list.M[["oo.o"]]), names))))
		xx <- xx[c(1,length(xx))]
		
		mmm <- matrix(c(
			diff(meancor(list.M[[paste0(nm1, ".", nm2, ".mp")]], names(i1), names(i2), nm=as.character(xx))),
			diff(meancor(list.M[[paste0(nm1, ".m")]], names(i1), names(i2), nm=as.character(xx))),
			diff(meancor(list.M[[paste0(nm1, ".", nm2, ".mm")]], names(i1), names(i2), nm=as.character(xx))),
			diff(meancor(list.M[[paste0(nm2, ".p")]], names(i1), names(i2), nm=as.character(xx))),
			diff(meancor(list.M[["oo.o"]], names(i1), names(i2), nm=as.character(xx))),
			diff(meancor(list.M[[paste0(nm2, ".m")]], names(i1), names(i2), nm=as.character(xx))),
			diff(meancor(list.M[[paste0(nm1, ".", nm2, ".pp")]], names(i1), names(i2), nm=as.character(xx))),
			diff(meancor(list.M[[paste0(nm1, ".p")]], names(i1), names(i2), nm=as.character(xx))),
			diff(meancor(list.M[[paste0(nm1, ".", nm2, ".pm")]], names(i1), names(i2), nm=as.character(xx)))
		), ncol=3)
		
		image(mmm[,ncol(mmm):1], asp=1, xaxt="n", yaxt="n", xlab="", ylab="", zlim=c(-0.5, 0.5), col=colorRampPalette(c("blue","white","red"))(1024))
		if (nm1 == "ie") {
			mtext(as.expression(phen.expression[names(i2)]), 2, line=1.5, xpd=NA, col=default.cols[names(i2)], cex=0.9)
		}
		if (nm2 == "st") {
			mtext(as.expression(phen.expression[names(i1)]), 1, line=1.5, xpd=NA, col=default.cols[names(i1)], cex=0.9)
		}
	}
	
	col.x <- col.sim[4]
	col.y <- col.sim[2]
	evol.dist <- 1.5
	plot(NULL, xlim=c(-2.5, 2.5), ylim=c(-2.5, 2.5), xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
	points(0, 0, pch=pch.sim["o"])
	points(x=c(0,0), y=evol.dist*c(-1,1), pch=pch.sim[c("m","p")], col=col.y, cex=2)
	arrows(x0=c(0,0), y0=c(0,0), x1=c(0,0), y1=0.8*evol.dist*c(-1,1), col=col.y, length=0.05)
	points(x=evol.dist*c(-1,1), y=c(0,0), pch=pch.sim[c("m","p")], col=col.x, cex=2)
	arrows(x0=c(0,0), y0=c(0,0), x1=0.8*evol.dist*c(-1,1), y1=c(0,0), col=col.x, length=0.05)
	points(x=evol.dist*c(-1,-1,1,1), y=evol.dist*c(-1,1,1,-1), pch=pch.sim[c("mm","mp","pp","pm")], col="darkgray", cex=2)
	arrows(x0=c(0,0,0,0), y0=c(0,0,0,0), x1=0.8*evol.dist*c(-1,-1,1,1), y1=0.8*evol.dist*c(-1,1,1,-1), col="darkgray",length=0.05)
	arrows(x0=-1.8, x1=2.2, y0=-2, length=0.05)
	text(0, -2.5, "Rob. comp. 1", xpd=NA)
	arrows(x0=-2.2, y0=-1.8, y1=1.8, length=0.05)
	text(-2.7, -0.2, "Rob. comp. 2", srt=90, xpd=NA)
	
	par(mar=c(1,5,1,1))
	image(x=0:1, y=seq(-0.5, 0.5, length.out=1025), z=t(matrix(seq(-0.5, 0.5, length.out=1024))), xlab="", ylab=expression(Delta*r), xaxt="n", zlim=c(-0.5, 0.5), col=colorRampPalette(c("blue","white","red"))(1024))

dev.off()
