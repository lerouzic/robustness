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
W0                <- matrix(rnorm(n.genes^2, sd=default.initsd), ncol=n.genes)
reps              <- default.sim.reps
test.rep          <- default.rob.reps
grad.effect       <- 0.01
N                 <- default.N 
G                 <- 10000
every             <- round(G/100)
force.run         <- !use.cache

avgcol <- function(c1, c2) rgb(colorRamp(c(c1,c2))(0.5), max=255)

col.sim <- c(oo="black", ie=COL.ENVCAN, le=COL.HOMEO, im=COL.GENCAN, lm=COL.SOM, st=COL.STAB, ie.lm=avgcol(COL.ENVCAN,COL.SOM), le.lm=avgcol(COL.HOMEO,COL.SOM), ie.st=avgcol(COL.ENVCAN,COL.STAB))
col.phen <- c(fitness="black", initenv=COL.ENVCAN, lateenv=COL.HOMEO, initmut=COL.GENCAN, latemut=COL.SOM, stability=COL.STAB)
lty.sim <- c(p=0, m=0, o=0, pp=0, pm=0, mp=0, mm=0)
pch.sim <- c(p=2, m=6, o=1, pp=24, mm=25, pm=0, mp=5)
max.points <- 20

phen <- c(
	fitness="Fitness",
    initenv=substitute(x~(y), list(x=TERM.ENVCAN.LONG, y=ABBRV.ENVCAN[[1]])),
    lateenv=substitute(x~(y), list(x=TERM.HOMEO.LONG, y=ABBRV.HOMEO[[1]])),
    initmut=substitute(x~(y), list(x=TERM.GENCAN.LONG, y=ABBRV.GENCAN[[1]])),
    latemut=substitute(x~(y), list(x=TERM.SOM.LONG, y=ABBRV.SOM[[1]])),
    stability=substitute(x~(y), list(x=TERM.STAB.LONG, y=ABBRV.STAB[[1]])))

allplots <- function(list.sim, what="fitness", xlim=NULL, ylim=NULL, xlab="Generations", ylab=what, lwd=3, FUN=mean, ...) {
	if(is.null(xlim)) xlim <- c(0, as.numeric(rev(names(list.sim[[1]]$mean))[1]))
	allv <- do.call(rbind, lapply(list.sim, function(x) do.call(rbind, lapply(x$mean, function(xx) FUN(xx[[what]])))))
	if(is.null(ylim)) ylim <- c(min(allv), max(allv))
	plot(NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab="", main=as.expression(phen[what]), col.main=col.phen[what], ...)
	for (nss in names(list.sim)) {
		nnss <- strsplit(nss, split="\\.")[[1]]
		if (length(nnss) == 3) nnss <- c(paste(nnss[1], nnss[2], sep="."), nnss[3])
		xpp <- seq_along(as.numeric(names(list.sim[[nss]]$mean)))
		xpp <- round(seq(xpp[1], xpp[length(xpp)], length.out=min(max.points, length(xpp))))
		lines(as.numeric(names(list.sim[[nss]]$mean))[xpp], sapply(list.sim[[nss]]$mean, function(x) FUN(x[[what]]))[xpp], col=col.sim[nnss[1]], bg=col.sim[nnss[1]], lty=lty.sim[nnss[2]], pch=pch.sim[nnss[2]], type="o", lwd=lwd)
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
# Some pairs of indexes:
# P1: low correlation, early environmental/late genetic
# P2: medium correlation, late environmental/late genetic
# P3: strong correlation, early environmental/stability

torun <- list(
	oo.o = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,0,0,0,0)), 
		reps=reps, series.name="figG-null", force.run=force.run), 
	ie.m = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(-grad.effect,0,0,0,0)), 
		reps=reps, series.name="figG-ie-m", force.run=force.run),
	le.m = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,-grad.effect,0,0,0)), 
		reps=reps, series.name="figG-le-m", force.run=force.run),
	im.m = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,0,0,-grad.effect,0)), 
		reps=reps, series.name="figG-lm-m", force.run=force.run),
	lm.m = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,0,0,-grad.effect,0)), 
		reps=reps, series.name="figG-lm-m", force.run=force.run),
	st.m = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,0,0,0,-grad.effect)), 
		reps=reps, series.name="figG-st-m", force.run=force.run),
	ie.p = function()pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(grad.effect,0,0,0,0)), 
		reps=reps, series.name="figG-ie-p", force.run=force.run),
	le.p = function()pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,grad.effect,0,0,0)), 
		reps=reps, series.name="figG-le-p", force.run=force.run),
	im.p = function()pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,0,grad.effect,0,0)), 
		reps=reps, series.name="figG-im-p", force.run=force.run),
	lm.p = function()pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,0,0,grad.effect,0)), 
		reps=reps, series.name="figG-lm-p", force.run=force.run),
	st.p = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,0,0,0,grad.effect)), 
		reps=reps, series.name="figG-st-p", force.run=force.run),
	ie.lm.mm = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(-grad.effect,0,0,-grad.effect,0)), 
		reps=reps, series.name="figI-ie-lm-mm", force.run=force.run),
	ie.lm.mp = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(-grad.effect,0,0,grad.effect,0)), 
		reps=reps, series.name="figI-ie-lm-mp", force.run=force.run),
	ie.lm.pm = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(grad.effect,0,0,-grad.effect,0)), 
		reps=reps, series.name="figI-ie-lm-pm", force.run=force.run),
	ie.lm.pp = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(grad.effect,0,0,grad.effect,0)), 
		reps=reps, series.name="figI-ie-lm-pp", force.run=force.run),
	le.lm.mm = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,-grad.effect,0,-grad.effect,0)), 
		reps=reps, series.name="figI-le-lm-mm", force.run=force.run),
	le.lm.mp = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,-grad.effect,0,grad.effect,0)), 
		reps=reps, series.name="figI-le-lm-mp", force.run=force.run),
	le.lm.pm = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,grad.effect,0,-grad.effect,0)), 
		reps=reps, series.name="figI-le-lm-pm", force.run=force.run),
	le.lm.pp = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(0,grad.effect,0,grad.effect,0)), 
		reps=reps, series.name="figI-le-lm-pp", force.run=force.run),
	ie.st.mm = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(-grad.effect,0,0,0,-grad.effect)), 
		reps=reps, series.name="figI-ie-st-mm", force.run=force.run),
	ie.st.mp = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(-grad.effect,0,0,0,grad.effect)), 
		reps=reps, series.name="figI-ie-st-mp", force.run=force.run),
	ie.st.pm = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(grad.effect,0,0,0,-grad.effect)), 
		reps=reps, series.name="figI-ie-st-pm", force.run=force.run),
	ie.st.pp = function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=every, grad.rob=c(grad.effect,0,0,0,grad.effect)), 
		reps=reps, series.name="figI-ie-st-pp", force.run=force.run)
)

list.sim <- mclapply(torun, function(ff) ff(), mc.cores=min(length(torun), ceiling(mc.cores/reps)))


pdf("figI.pdf", width=8, height=12)
	layout(rbind(1:2,3:4,5:6))
	par(cex=1)
	
	allboxes(list.sim, "ie","lm",what="initenv", G="5000", lwd=3)
	allboxes(list.sim, "ie","lm",what="latemut", G="5000", lwd=3)
	
	allboxes(list.sim, "le","lm",what="lateenv", G="5000", lwd=3)
	allboxes(list.sim, "le","lm",what="latemut", G="5000", lwd=3)
	
	allboxes(list.sim, "ie","st",what="initenv", G="5000", lwd=3)
	allboxes(list.sim, "ie","st",what="stability", G="5000", lwd=3)
	
dev.off()
	
