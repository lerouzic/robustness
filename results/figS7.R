#!/usr/bin/env Rscript

source("./commonsim.R")
source("./terminology.R")
source("./defaults.R")

################## Options

n.genes          <- default.n
sel.genes        <- default.nsel
s                <- c(rep(default.s, sel.genes), rep(0, n.genes-sel.genes))
W0               <- NA
reps             <- default.sim.reps
test.rep         <- default.rob.reps
N                <- default.N
G                <- 20000
mut.rate         <- default.mut.rate
force.run        <- !use.cache

mc.cores         <- default.mc.cores

defaults         <- c(m=mut.rate, N=N, g=n.genes, sg=sel.genes, s=default.s)

nb.values        <- 11
mut.values       <- 10^seq(-3,-1, length.out=nb.values)
N.values         <- round(10^seq(1, 4, length.out=nb.values))
genes.values     <- round(seq(3, 20, length.out=nb.values))
selg.values      <- 1:6
s.values         <- 10^seq(-1, 2, length.out=nb.values)

phen             <- c(list(fitness="Fitness"), phen.expression)

ylims            <- list(
                       fitness=c(0.85,1), 
                       initenv=c(-35,-15), 
                       lateenv=c(-15,0), 
                       initmut=c(-13,-4), 
                       latemut=c(-15,-5), 
                       stability=c(-35,-12))

captions <- c(
	m=expression("Mutation rate ("*mu*")"),
	N="Population size (N)",
	g="Number of genes (n)",
	sg="Selected genes (n\')",
	s="Selection strength (s)")


###################### Functions
torun.mut <- lapply(mut.values, function(mm) 
	substitute(function() sim.run.reps(
		W0, 
		list(s=s, G=G, N=N, rep=test.rep, summary.every=G, mut.rate=mm), 
		reps=reps, series.name=paste0("figL-ref-m", signif(mm, digits=2)), 
		force.run=force.run, mc.cores=min(reps, mc.cores)), 
	list(mm=mm)))
	
torun.N <- lapply(N.values, function(nn) 
	substitute(function() sim.run.reps(
		W0, 
		list(s=s, G=G, N=nn, rep=test.rep, summary.every=G, mut.rate=mut.rate), 
		reps=reps, series.name=paste0("figL-ref-N", nn), 
		force.run=force.run, mc.cores=min(reps, mc.cores)), 
	list(nn=nn)))
	
torun.genes <- lapply(genes.values, function(gg) 
	substitute(function() sim.run.reps(
		W0, 
		list(s=c(rep(10, min(gg, sel.genes)), rep(0 ,gg-min(gg,sel.genes))), theta=rep(NA, gg), G=G, N=N, rep=test.rep, summary.every=G, mut.rate=mut.rate), 
		reps=reps, series.name=paste0("figL-ref-g", gg), 
		force.run=force.run, mc.cores=min(reps, mc.cores)), 
	list(gg=gg)))
	
torun.selg <- lapply(selg.values, function(sg)
	substitute(function() sim.run.reps(
		W0, 
		list(s=c(rep(10,sg),rep(0,n.genes-sg)), G=G, N=N, rep=test.rep, summary.every=G, mut.rate=mut.rate), 
		reps=reps, series.name=paste0("figL-ref-sg", sg), 
		force.run=force.run, mc.cores=min(reps, mc.cores)), 
	list(sg=sg)))
	
torun.sel <- lapply(s.values, function(ss)
	substitute(function() sim.run.reps(
		W0, 
		list(s=c(rep(ss,sel.genes),rep(0,n.genes-sel.genes)), G=G, N=N, rep=test.rep, summary.every=G, mut.rate=mut.rate), 
		reps=reps, series.name=paste0("figL-ref-s", signif(ss, digits=2)), 
		force.run=force.run, mc.cores=min(reps, mc.cores)), 
	list(ss=ss)))

torun <- setNames(
	c(torun.mut, torun.N, torun.genes, torun.selg, torun.sel),
	c(	paste0("ref.m", signif(mut.values, digits=2)), 
		paste0("ref.N", N.values), 
		paste0("ref.g", genes.values), 
		paste0("ref.sg", selg.values), 
		paste0("ref.s", signif(s.values, digits=2))))
	
	
################# Calc	
list.sim <- mclapply(torun, function(ff) eval(ff)(), mc.cores=min(length(torun), ceiling(mc.cores/reps)))

################# Figure

ww <- c("fitness", "initenv", "lateenv","initmut","latemut", "stability")

pdf("figS7.pdf", width=10, height=14)
	layout(matrix(1:(length(ww)*length(captions)), ncol=length(captions)))
	par(mar=c(0.5, 0.5, 0.5, 0.1), oma=c(5, 4, 0, 0))
	for (pp in names(captions)) {
		for (what in ww) {
			ls <- list.sim[grep(names(list.sim), pattern=paste0("ref\\.",pp,"\\d"))]
			xval <- sapply(strsplit(names(ls), split=paste0("\\.",pp)), function(sp) as.numeric(sp[2]))
			yval <- sapply(ls, function(x) mean(x$mean[[as.character(G)]][[what]]))
			yvar <- sapply(ls, function(x) mean(x$var[[as.character(G)]][[what]]))
			plot(NULL, log=if(pp %in% c("g","sg")) "" else "x", xaxt="n", yaxt="n", xlab="", ylab="", xlim=range(xval), ylim=ylims[[what]])
			arrows(x0=xval, y0=yval-sqrt(yvar), y1=yval+sqrt(yvar), code=3, length=0, angle=90, col="darkgray")
			points(xval, yval, type="p")
			
			abline(v=defaults[pp], col="gray", lty=3, lwd=3)
			if(pp==names(captions)[1]) {
				mtext(2, text=as.expression(phen[[what]]), line=2, xpd=NA)
				axis(2)
			}
			if (what==ww[length(ww)]) {
				mtext(1, text=captions[pp], line=3, xpd=NA)
				axis(1)
			}
		}
	}
dev.off()
