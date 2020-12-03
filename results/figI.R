#!/usr/bin/env Rscript

source("./commonpure.R")
source("./terminology.R")
source("./defaults.R")

library(parallel)
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
G                 <- 10000
every             <- round(G/100)
force.run         <- !use.cache

avgcol <- function(c1, c2) rgb(colorRamp(c(c1,c2))(0.5), max=255)

phen <- c(list(fitness="Fitness"), phen.expression)
col.sim <- c(oo="black", ie=COL.ENVCAN, le=COL.HOMEO, im=COL.GENCAN, lm=COL.SOM, st=COL.STAB, ie.lm=avgcol(COL.ENVCAN,COL.SOM), le.lm=avgcol(COL.HOMEO,COL.SOM), ie.st=avgcol(COL.ENVCAN,COL.STAB), ie.im=avgcol(COL.ENVCAN, COL.GENCAN))
col.phen <- c(fitness="black", initenv=COL.ENVCAN, lateenv=COL.HOMEO, initmut=COL.GENCAN, latemut=COL.SOM, stability=COL.STAB)
lty.sim <- c(p=0, m=0, o=0, pp=0, pm=0, mp=0, mm=0)
pch.sim <- c(p=2, m=6, o=1, pp=24, mm=25, pm=0, mp=5)
max.points <- 20

allplots <- function(list.sim, r1, r2, what="fitnesses", ylab=as.expression(phen[what]), r1.series="p", r2.series="m", comb="pm", relative=NULL, xlim=NULL, ylim=NULL, xlab="Generations", lwd=1, ...) {
	G <- as.numeric(names(list.sim[[1]]$mean))
	which.G <- seq(1, length(G), length.out=20)
	mylist.sim <- lapply(list.sim, function(x) sapply(x$mean, function(xx) mean(xx[[what]])))
	mylist.sim <- mylist.sim[c("oo.o", paste(r1, r1.series,sep="."), paste(r2, r2.series,sep="."), paste(r1, r2, comb, sep="."))]
	if (length(relative) > 0) {
		mylist.sim <- lapply(mylist.sim, function(xx) xx-mean(mylist.sim[[relative]]))
		mylist.sim[[relative]] <- NULL
	}
	if (is.null(xlim)) xlim <- range(G)
	if (is.null(ylim)) ylim <- range(unlist(mylist.sim), na.rm=TRUE)
	plot(NULL, xlim=xlim, ylim=ylim, ylab=ylab, xlab=xlab, ...)
	for (ni in names(mylist.sim)) {
		sp <- strsplit(ni, split="\\.")[[1]]
		sp1 <- if(length(sp) == 2) sp[1] else paste(sp[1],sp[2],sep=".")
		sp2 <- sp[length(sp)]
		points(G[which.G], mylist.sim[[ni]][which.G], pch=pch.sim[sp2], col=col.sim[sp1], lwd=lwd)
	}
}

allboxes <- function(list.sim, r1, r2, G=NULL, what="fitnesses", xlab=phen[what], relative=NULL, ...) {
	if (is.null(G)) G <- rev(names(list.sim[[1]]$mean[[1]]))[1]
	mylist.sim <- lapply(list.sim, function(x) sapply(x$full, function(xx) mean(xx[[G]][[what]])))
	mylist.sim <- mylist.sim[c("oo.o", paste(r1, c("m","p"),sep="."), paste(r2, c("m","p"),sep="."), paste(r1, r2, c("mm","mp","pm","pp"), sep="."))]
	if (length(relative) > 0) {
		mylist.sim <- lapply(mylist.sim, function(xx) xx-mean(mylist.sim[[relative]]))
		mylist.sim[[relative]] <- NULL
	}
	boxplot(mylist.sim, horizontal=TRUE, xlab=as.expression(xlab), at=-c(if(length(relative) == 0) 1 else NULL, 3:4, 6:7, 9:12), yaxt="n", border=col.sim[sapply(strsplit(names(mylist.sim), split="\\."), function(ss) if(length(ss) == 2) ss[1] else paste(ss[1],ss[2],sep="."))], frame=FALSE,  ...)
	labels <- if (col.phen[what] == col.sim[r1])
			c(if(length(relative) == 0) "Control" else NULL, "Direct -", "Direct +", "Indirect -", "Indirect +", "D-; I-", "D-; I+", "D+; I-", "D+; I+")
		else
			c(if(length(relative) == 0) "Control" else NULL, "Indirect -", "Indirect +", "Direct -", "Direct +", "D-; I-", "D+; I-", "D-; I+", "D+; I+")
	axis(2, las=2, at=-c(if(length(relative) == 0) 1 else NULL, 3:4, 6:7, 9:12), tick=FALSE, lty=0, line=-1, label=labels)
}

geom.analysis <- function(list.sim, r1, r2, what1, what2, G=NULL, xlim=NULL, ylim=NULL, plot.grad=TRUE, means=TRUE, ..., ccol=setNames(rainbow(8), c("p1","pm","m2","mm","m1","mp","p2","pp"))) {
	if (is.null(G)) G <- rev(names(list.sim[[1]]$mean))[1]
	min.G <- names(list.sim[[1]]$mean)[1]
				
	O.x <- mean(list.sim[["oo.o"]]$mean[[min.G]][[what1]])
	O.y <- mean(list.sim[["oo.o"]]$mean[[min.G]][[what2]])
	
	ssim <- c(paste(r1, c("m","p"),sep="."), paste(r2, c("m","p"),sep="."), paste(r1, r2, c("mm","mp","pm","pp"), sep="."))
	mylist.sim1 <- lapply(list.sim[ssim], function(x) sapply(x$full, function(xx) mean(xx[[G]][[what1]])))
	mylist.sim2 <- lapply(list.sim[ssim], function(x) sapply(x$full, function(xx) mean(xx[[G]][[what2]])))
	
	names(mylist.sim1) <- names(mylist.sim2) <- c("m1","p1", "m2","p2","mm", "mp", "pm", "pp")
	
	if (means) {
		mylist.sim1 <- lapply(mylist.sim1, mean)
		mylist.sim2 <- lapply(mylist.sim2, mean)
	}
	
	xlim <- if (is.null(xlim)) range(unlist(mylist.sim1))
	if (is.null(ylim)) ylim <- range(unlist(mylist.sim2))
	

	
	plot(NULL, xlim=xlim, ylim=ylim, xlab=as.expression(phen[what1]), ylab=as.expression(phen[what2]), ...)
	points(O.x, O.y, pch=1, cex=3, col="black")
	if (plot.grad) {
		grad.size.x <- abs(min(O.x-xlim[1], xlim[2]-O.x))
		grad.size.y <- abs(min(O.y - ylim[1], ylim[2]-O.y))
		arrows(x0=O.x, y0=O.y, x1=O.x+c(1, 1, 0, -1, -1, -1, 0, 1)*grad.size.x, y1=O.y+c(0, -1, -1, -1, 0, 1, 1, 1, 0)*grad.size.y, col=ccol)
	}
	points(unlist(mylist.sim1), unlist(mylist.sim2), col=ccol[rep(names(mylist.sim1), sapply(mylist.sim1, length))])
	text(sapply(mylist.sim1, mean), sapply(mylist.sim2, mean), names(mylist.sim1), pos=ifelse(sapply(mylist.sim1, mean) > O.x, 2, 4))
}


rel.response.onesim <- function(mm, mp, pm, pp, what1, what2, G=NULL, normalize=FALSE) {
	if (is.null(G)) G <- rev(names(mm))[1]
	
	allpoints <- rbind( mm = c(mean(mm[[G]][[what1]]), mean(mm[[G]][[what2]])), 
						mp = c(mean(mp[[G]][[what1]]), mean(mp[[G]][[what2]])), 
						pm = c(mean(pm[[G]][[what1]]), mean(pm[[G]][[what2]])), 
						pp = c(mean(pp[[G]][[what1]]), mean(pp[[G]][[what2]])))
						
	if (normalize) {
		allpoints[,1] <- allpoints[,1]-allpoints["mm",1]
		allpoints[,2] <- allpoints[,2]-allpoints["mm",2]
		allpoints[,1] <- allpoints[,1] / allpoints["pp",1]
		allpoints[,2] <- allpoints[,2] / allpoints["pp",2]
	}
	
	# Direction of the most evolvability: mm to pp
	a <- (allpoints["pp",2] - allpoints["mm",2])/(allpoints["pp",1] - allpoints["mm",1])
	ap <- - 1/a # perpendicular direction
	
	mpp.x <- (allpoints["mp",1] + ap*allpoints["mp",2])/(1+ap^2)
	mpp.y <- ap*(allpoints["mp",1] + ap*allpoints["mp",2])/(1+ap^2)
	
	pmp.x <- (allpoints["pm",1] + ap*allpoints["pm",2])/(1+ap^2)
	pmp.y <- ap*(allpoints["pm",1] + ap*allpoints["pm",2])/(1+ap^2)
	
	dist.mm.pp <- sqrt((allpoints["mm",1]-allpoints["pp",1])^2 + (allpoints["mm",2] - allpoints["pp",2])^2)
	dist.mp.pm <- sqrt((mpp.x - pmp.x)^2 + (mpp.y - pmp.y)^2)
	return(c(max.dist=dist.mm.pp, min.dst=dist.mp.pm, ratio=dist.mp.pm/dist.mm.pp))		 
}

rel.response <- function(list.sim, r1, r2, what1, what2, G=NULL, normalize=FALSE) {
	
	ssim <- c(paste(r1, r2, c("mm","mp","pm","pp"), sep="."))
		
	if (any(!ssim %in% names(list.sim))) return(list(mean.ratio=NA, sd.ratio=NA, se.ratio=NA))
	
	# loop over replicates
	rep <- length(list.sim[[ssim[1]]]$full)
	ans <- lapply(1:rep, function(i) 
		rel.response.onesim(mm = list.sim[[ssim[1]]]$full[[i]], 
		                    mp = list.sim[[ssim[2]]]$full[[i]], 
		                    pm = list.sim[[ssim[3]]]$full[[i]], 
		                    pp = list.sim[[ssim[4]]]$full[[i]],
		                    what1=what1, what2=what2, G=G, normalize=normalize))
	
	rat <- sapply(ans, "[", "ratio.mp")
	list(mean.ratio = mean(rat), sd.ratio = sd(rat), se.ratio = sd(rat)/sqrt(rep))
}

indx <- c(ie=1, le=2, im=3, lm=4, st=5)
gradvec1 <- function(i, grd) c(rep(0, i-1), grd, rep(0, length(indx)-i))
gradvec2 <- function(i1, i2, grd1, grd2) c(rep(0, i1-1), grd1, rep(0, i2-i1-1), grd2, rep(0, length(indx)-i2))


# Make a list of simulation names and gradients

torun <- list(
	oo.o = list(series.name="figG-null", grad.rob=c(0,0,0,0,0)))
for (i in indx) {
	nm <- names(indx)[i]
	torun[[paste(nm, "m", sep=".")]] <- list(series.name=paste("figG", nm, "m", sep="-"), grad.rob=gradvec1(i, -grad.effect))
	torun[[paste(nm, "p", sep=".")]] <- list(series.name=paste("figG", nm, "p", sep="-"), grad.rob=gradvec1(i,  grad.effect))
}
for (i1 in 1:(length(indx)-1))
	for (i2 in (i1+1):length(indx)) {
		nm1 <- names(indx)[i1]
		nm2 <- names(indx)[i2]
		torun[[paste(nm1, nm2, "mm", sep=".")]] <- list(series.name=paste("figI", nm1, nm2, "mm", sep="-"), grad.rob=gradvec2(i1, i2, -grad.effect, -grad.effect))
		torun[[paste(nm1, nm2, "mp", sep=".")]] <- list(series.name=paste("figI", nm1, nm2, "mp", sep="-"), grad.rob=gradvec2(i1, i2, -grad.effect,  grad.effect))
		torun[[paste(nm1, nm2, "pm", sep=".")]] <- list(series.name=paste("figI", nm1, nm2, "pm", sep="-"), grad.rob=gradvec2(i1, i2,  grad.effect, -grad.effect))
		torun[[paste(nm1, nm2, "pp", sep=".")]] <- list(series.name=paste("figI", nm1, nm2, "pp", sep="-"), grad.rob=gradvec2(i1, i2,  grad.effect,  grad.effect))
	}

list.sim <- mclapply(torun, function(ff) 
	pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=ff$grad.rob),	reps=reps, series.name=ff$series.name, force.run=force.run), 
	mc.cores=min(length(torun), ceiling(mc.cores/reps)))


indwhat <- c(ie="initenv", le="lateenv", im="initmut", lm="latemut", st="stability")
tt <- matrix("", ncol=length(indwhat)-1, nrow=length(indwhat)-1)
colnames(tt) <- names(indwhat)[1:(length(indwhat)-1)]
rownames(tt) <- names(indwhat)[2:length(indwhat)]

for (i1 in 1:(length(indwhat)-1))
	for (i2 in(i1+1):length(indwhat)) {
		n1 <- names(indwhat)[i1]
		n2 <- names(indwhat)[i2]
		rr <- rel.response(list.sim, n1, n2, indwhat[i1], indwhat[i2], normalize=FALSE)
		tt[n2, n1] <- paste0(round(rr$mean.ratio, digits=2), " +/- ", round(rr$sd.ratio, digits=2))
	}



geom.analysis(list.sim, "le","lm", what1="lateenv", what2="latemut")
rel.response(list.sim, "le","lm", what1="lateenv", what2="latemut")



pdf("figI.pdf", width=8, height=12) 
	layout(rbind(1:2,3:4,5:6,7:8))
	par(cex=1)
	
	allplots(list.sim, "ie","lm", what="initenv", r1.series="p", r2.series="m", comb="pm", lwd=3)
	allplots(list.sim, "ie","lm", what="latemut", r1.series="p", r2.series="m", comb="pm", lwd=3)
	
	allplots(list.sim, "le","lm", what="lateenv", r1.series="p", r2.series="m", comb="pm", lwd=3)
	allplots(list.sim, "le","lm", what="latemut", r1.series="p", r2.series="m", comb="pm", lwd=3)
	
#~ 	allplots(list.sim, "ie","st", what="initenv", r1.series="p", r2.series="m", comb="pm", lwd=3)
#~ 	allplots(list.sim, "ie","st", what="stability", r1.series="p", r2.series="m", comb="pm", lwd=3)

	allplots(list.sim, "ie","st", what="initenv", r1.series="m", r2.series="p", comb="mp", lwd=3)
	allplots(list.sim, "ie","st", what="stability", r1.series="m", r2.series="p", comb="mp", lwd=3)

	allplots(list.sim, "ie","im", what="initenv", r1.series="m", r2.series="p", comb="mp", lwd=3)
	allplots(list.sim, "ie","im", what="initmut", r1.series="m", r2.series="p", comb="mp", lwd=3)
dev.off()


pdf("figIb.pdf", width=8, height=12)
	layout(rbind(1:2,3:4,5:6))
	par(cex=1)
	
	allboxes(list.sim, "ie","lm",what="initenv", G="5000", lwd=3, ylim=c(-40, -5))
	allboxes(list.sim, "ie","lm",what="latemut", G="5000", lwd=3, ylim=c(-12, -6))
	
	allboxes(list.sim, "le","lm",what="lateenv", G="5000", lwd=3, ylim=c(-15, -4))
	allboxes(list.sim, "le","lm",what="latemut", G="5000", lwd=3, ylim=c(-15, -6))
	
	allboxes(list.sim, "ie","st",what="initenv", G="5000", lwd=3)
	allboxes(list.sim, "ie","st",what="stability", G="5000", lwd=3)
	
dev.off()
	

# Explanatory figure for the relative response in the direction of most evolutionary
# resistance
# Note: this scheme is only partially symbolic, mm and pp are expected to be centered around (0,0)

library(ellipse)
pdf("figIe.pdf", width=5, height=5)
par(mar=c(2,2,0.1,0.1))
plot(NULL, xlim=c(-1,1), ylim=c(-1,1), xaxt="n", yaxt="n", xlab="", ylab="",asp=1)

mtext(expression(rho[1]), 1, line=1)
mtext(expression(rho[2]), 2, line=1)

pp <- c(0.8,0.8)
mm <- c(-0.8,-0.8)
pm <- c(0.4, -0.1)
mp <- c(-0.45, 0.)
points(c(pp[1], mm[1]), c(pp[2], mm[2]), cex=2, pch=19, col=c("blue","red"))
text(c(pp[1], mm[1]), c(pp[2], mm[2]), c(expression(rho[1]*"+,"*rho[2]*"+"), expression(rho[1]*"-,"*rho[2]*"-")), col=c("blue","red"), pos=c(2,4))

slp <- (pp[2]-mm[2])/(pp[1]-mm[1])
abline(a=0, b=slp, lty=3, col="darkgray")
abline(a=0, b=-1/slp, lty=2, col="darkgray")
sq.sz <- 0.05
arrows(	x0=c(-slp*sq.sz, (slp-1)*sq.sz)/sqrt(1+slp^2),
		y0=c(sq.sz, (1+slp)*sq.sz)/sqrt(1+slp^2),
		x1=c((slp-1)*sq.sz, sq.sz)/sqrt(1+slp^2), 
		y1=c((1+slp)*sq.sz, slp*sq.sz)/sqrt(1+slp^2), 
		lty=1, col="darkgray", code=0)

points(c(pm[1], mp[1]), c(pm[2], mp[2]), cex=2, pch=19, col=c("darkgreen", "orange"))
text(c(pm[1], mp[1]), c(pm[2], mp[2]), c(expression(rho[1]*"+,"*rho[2]*"-"), expression(rho[1]*"-,"*rho[2]*"+")), col=c("darkgreen","orange"), pos=c(4,2))


pm.proj <- c((pm[1]-pm[2]/slp)/(1+1/slp^2), -(pm[1]-pm[2]/slp)/slp/(1+1/slp^2))
mp.proj <- c((mp[1]-mp[2]/slp)/(1+1/slp^2), -(mp[1]-mp[2]/slp)/slp/(1+1/slp^2))

arrows(x0=c(pm[1], mp[1]), y0=c(pm[2], mp[2]), x1=c(pm.proj[1], mp.proj[1]) , y1=c(pm.proj[2], mp.proj[2]), col="red", length=0.1)

arrows(x0=pp[1], x1=mm[1], y0=pp[2], y1=mm[2], code=3, length=0.2, angle=90, col="blue")
text(-0.5*slp, -0.5*slp, expression(Delta[1]), col="blue", pos=1)

arrows(x0=mp.proj[1], x1=pm.proj[1], y0=mp.proj[2], y1=pm.proj[2], col="deeppink", code=3, length=0.2, angle=90)
text(-0.2/slp, 0.2/slp, expression(Delta[2]), col="deeppink", pos=4)

text(0.3, -0.7, expression(R==Delta[2]/Delta[1]))

lines(ellipse(0.8, scale=0.4*c(1,1)), col="black", lty=3)

dev.off()
