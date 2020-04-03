#!/usr/bin/env Rscript

source("./commonpure.R")
source("./terminology.R")
source("./defaults.R")

library(parallel)
mc.cores <- default.mc.cores

use.cache <- TRUE

n.genes          <- default.n
sel.genes        <- default.nsel
s                <- c(rep(default.s, sel.genes), rep(0, n.genes-sel.genes))
W0               <- matrix(rnorm(n.genes^2, sd=default.initsd), ncol=n.genes)
reps             <- default.sim.reps
test.rep         <- default.rob.reps
N                <- default.N
G                <- 20000
mut.rate         <- default.mut.rate
force.run        <- !use.cache

nb.values <- 11

defaults <- c(m=mut.rate, N=N, g=n.genes, sg=sel.genes, s=default.s)

mut.values <- 10^seq(-4,-1, length.out=nb.values)
N.values   <- round(10^seq(1, 4, length.out=nb.values))
genes.values <- round(seq(3, 20, length.out=nb.values))
selg.values  <- 1:6
s.values <- 10^seq(-2, 3, length.out=nb.values)

phen <- c(
	fitness="Fitness",
    initenv=substitute(x~(y), list(x=TERM.ENVCAN.SHORT, y=ABBRV.ENVCAN[[1]])),
    lateenv=substitute(x~(y), list(x=TERM.HOMEO.SHORT, y=ABBRV.HOMEO[[1]])),
    initmut=substitute(x~(y), list(x=TERM.GENCAN.SHORT, y=ABBRV.GENCAN[[1]])),
    latemut=substitute(x~(y), list(x=TERM.SOM.SHORT, y=ABBRV.SOM[[1]])),
    stability=substitute(x~(y), list(x=TERM.STAB.SHORT, y=ABBRV.STAB[[1]])))
    
captions <- c(
	m=expression("Mutation rate ("*mu*")"),
	N="Population size (N)",
	g="Number of genes (n)",
	sg="Selected genes (n\')",
	s="Selection strength (s)")

torun.mut <- lapply(mut.values, function(mm) 
	substitute(function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=G, mut.rate=mm), 
		reps=reps, series.name=paste0("figL-ref-m", signif(mm, digits=2)), force.run=force.run), list(mm=mm)))
torun.N <- lapply(N.values, function(nn) 
	substitute(function() pure.run.reps(W0, list(s=s, G=G, N=nn, rep=test.rep, summary.every=G, mut.rate=mut.rate), 
		reps=reps, series.name=paste0("figL-ref-N", nn), force.run=force.run), list(nn=nn)))
torun.genes <- lapply(genes.values, function(gg) 
	substitute(function() pure.run.reps(W0=matrix(rnorm(gg^2, sd=0.000001), ncol=gg), list(s=c(rep(10, sel.genes), rep(0 ,gg-sel.genes)), G=G, N=N, rep=test.rep, summary.every=G, mut.rate=mut.rate), 
		reps=reps, series.name=paste0("figL-ref-g", gg), force.run=force.run), list(gg=gg)))
torun.selg <- lapply(selg.values, function(sg)
	substitute(function() pure.run.reps(W0, list(s=c(rep(10,sg),rep(0,n.genes-sg)), G=G, N=N, rep=test.rep, summary.every=G, mut.rate=mut.rate), 
		reps=reps, series.name=paste0("figL-ref-sg", sg), force.run=force.run), list(sg=sg)))
torun.sel <- lapply(s.values, function(ss)
	substitute(function() pure.run.reps(W0, list(s=c(rep(ss,sel.genes),rep(0,n.genes-sel.genes)), G=G, N=N, rep=test.rep, summary.every=G, mut.rate=mut.rate), 
		reps=reps, series.name=paste0("figL-ref-s", signif(ss, digits=2)), force.run=force.run), list(ss=ss)))

torun <- setNames(
	c(torun.mut, torun.N, torun.genes, torun.selg, torun.sel),
	c(paste0("ref.m", signif(mut.values, digits=2)), paste0("ref.N", N.values), paste0("ref.g", genes.values), paste0("ref.sg", selg.values), paste0("ref.s", signif(s.values, digits=2))))
list.sim <- mclapply(torun, function(ff) eval(ff)(), mc.cores=min(length(torun), ceiling(mc.cores/reps)))

ylims <- list(initenv=c(-44,-15), lateenv=c(-30,0), initmut=c(-28,-4), latemut=c(-30,-5), stability=c(-44,-12))

pdf("figL.pdf", width=10, height=10)
	layout(matrix(1:25, ncol=5))
	par(mar=c(0.5, 0.5, 0.1, 0.1), oma=c(5, 4, 0, 0), xpd=NA)
	for (pp in names(captions)) {
		for (what in c("initenv", "lateenv","initmut","latemut", "stability")) {
			ls <- list.sim[grep(names(list.sim), pattern=paste0("ref\\.",pp,"\\d"))]
			xval <- sapply(strsplit(names(ls), split=paste0("\\.",pp)), function(sp) as.numeric(sp[2]))
			yval <- sapply(ls, function(x) mean(x$mean[[as.character(G)]][[what]]))
			yvar <- sapply(ls, function(x) mean(x$var[[as.character(G)]][[what]]))
			plot(NULL, log=if(pp %in% c("g","sg")) "" else "x", xaxt="n", yaxt="n", xlab="", ylab="", xlim=range(xval), ylim=ylims[[what]])
			arrows(x0=xval, y0=yval-sqrt(yvar), y1=yval+sqrt(yvar), code=3, length=0, angle=90, col="darkgray")
			points(xval, yval, type="o")
			
			abline(v=defaults[pp], col="gray", lty=3)
			if(pp==names(captions)[1]) {
				mtext(2, text=as.expression(phen[[what]]), line=2)
				axis(2)
			}
			if (what=="stability") {
				mtext(1, text=captions[pp], line=3)
				axis(1)
			}
		}
	}
dev.off()
