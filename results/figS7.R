#!/usr/bin/env Rscript

source("./commonsim.R")
source("./terminology.R")
source("./defaults.R")

################## Options

param <- default

param$s          <- c(rep(param$s, param$nsel), rep(0, param$n-param$nsel))
param$W0         <- NA
param$G          <- 5000
param$sim.reps   <- 20

cache.tag        <- "figL"

nb.values        <- 11
mut.values       <- 10^seq(-3,-1, length.out=nb.values)
smut.values      <- 10^seq(-3, 1, length.out=nb.values)
N.values         <- round(10^seq(1, 4, length.out=nb.values))
a.values         <- seq(0.1, 0.5, length.out=nb.values)
genes.values     <- round(seq(3, 20, length.out=nb.values))
selg.values      <- 1:6
d.values         <- seq(0.2, 1, length.out=nb.values)
s.values         <- 10^seq(-1, 2, length.out=nb.values)

phen             <- c(list(fitness="Fitness"), phen.expression)

ylims            <- list(
                       fitness=c(0.85,1), 
                       initenv=c(-35,-5), 
                       lateenv=c(-15,0), 
                       initmut=c(-17,0), 
                       latemut=c(-20,-0), 
                       stability=c(-35,-5))

captions <- c(
	m=expression("Mutation rate ("*nu*")"),
	sm=expression("Mutation size ("*sigma[nu]*")"),
	N="Population size (N)",
	a="Constitutive expression (a)",
	g="Number of genes (n)",
	sg="Selected genes (n\')",
	d="Network density (d)",
	s="Selection strength (s)")


###################### Listing simulations to run

torun.mut <- lapply(mut.values, function(mm) 
	substitute(function() sim.run.reps(
		param$W0, 
		list(s=param$s, G=param$G, N=param$N, rep=param$rob.reps, summary.every=param$G, mut.rate=mm), 
		reps=param$sim.reps, series.name=paste0(cache.tag, "-m", signif(mm, digits=2)), 
		force.run=!param$use.cache, mc.cores=min(param$sim.reps, param$mc.cores)), 
	list(mm=mm)))
	
torun.smut <- lapply(smut.values, function(sm) 
	substitute(function() sim.run.reps(
		param$W0, 
		list(s=param$s, G=param$G, N=param$N, rep=param$rob.reps, summary.every=param$G, mut.rate=param$mut.rate, mut.sd=sm), 
		reps=param$sim.reps, series.name=paste0(cache.tag, "-sm", signif(sm, digits=2)), 
		force.run=!param$use.cache, mc.cores=min(param$sim.reps, param$mc.cores)), 
	list(sm=sm)))
	
torun.a <- lapply(a.values, function(aa) 
	substitute(function() sim.run.reps(
		param$W0, 
		list(s=param$s, G=param$G, N=param$N, a=aa, rep=param$rob.reps, summary.every=param$G, mut.rate=param$mut.rate), 
		reps=param$sim.reps, series.name=paste0(cache.tag, "-a", signif(aa, digits=2)), 
		force.run=!param$use.cache, mc.cores=min(param$sim.reps, param$mc.cores)), 
	list(aa=aa)))
	
torun.N <- lapply(N.values, function(nn) 
	substitute(function() sim.run.reps(
		param$W0, 
		list(s=param$s, G=param$G, N=nn, rep=param$rob.reps, summary.every=param$G, mut.rate=param$mut.rate), 
		reps=param$sim.reps, series.name=paste0(cache.tag, "-N", nn), 
		force.run=!param$use.cache, mc.cores=min(param$sim.reps, param$mc.cores)), 
	list(nn=nn)))
	
torun.genes <- lapply(genes.values, function(gg) 
	substitute(function() sim.run.reps(
		param$W0, 
		list(
			s=c(rep(10, min(gg, param$nsel)), rep(0 ,gg-min(gg,param$nsel))), theta=rep(NA, gg), 
			G=param$G, N=param$N, rep=param$rob.reps, summary.every=param$G, mut.rate=param$mut.rate), 
		reps=param$sim.reps, series.name=paste0(cache.tag, "-g", gg), 
		force.run=!param$use.cache, mc.cores=min(param$sim.reps, param$mc.cores)), 
	list(gg=gg)))
	
torun.selg <- lapply(selg.values, function(sg)
	substitute(function() sim.run.reps(
		param$W0, 
		list(s=c(rep(10,sg),rep(0,param$n-sg)), G=param$G, N=param$N, rep=param$rob.reps, summary.every=param$G, mut.rate=param$mut.rate), 
		reps=param$sim.reps, series.name=paste0(cache.tag, "-sg", sg), 
		force.run=!param$use.cache, mc.cores=min(param$sim.reps, param$mc.cores)), 
	list(sg=sg)))
	
torun.d <- lapply(d.values, function(dd) 
	substitute(function() sim.run.reps(
		param$W0, 
		list(s=param$s, G=param$G, N=param$N, density=dd, rep=param$rob.reps, summary.every=param$G, mut.rate=param$mut.rate), 
		reps=param$sim.reps, series.name=paste0(cache.tag, "-d", signif(dd, digits=2)), 
		force.run=!param$use.cache, mc.cores=min(param$sim.reps, param$mc.cores)), 
	list(dd=dd)))
	
torun.sel <- lapply(s.values, function(ss)
	substitute(function() sim.run.reps(
		param$W0, 
		list(
			s=c(rep(ss,param$nsel),rep(0,param$n-param$nsel)), 
			G=param$G, N=param$N, rep=param$rob.reps, summary.every=param$G, mut.rate=param$mut.rate), 
		reps=param$sim.reps, series.name=paste0(cache.tag, "-s", signif(ss, digits=2)), 
		force.run=!param$use.cache, mc.cores=min(param$sim.reps, param$mc.cores)), 
	list(ss=ss)))

torun <- setNames(
	c(torun.mut, torun.smut, torun.a, torun.N, torun.genes, torun.selg, torun.d, torun.sel),
	c(	paste0("ref.m", signif(mut.values, digits=2)), 
		paste0("ref.sm", signif(smut.values, digits=2)), 
		paste0("ref.a", signif(a.values, digits=2)),
		paste0("ref.N", N.values), 
		paste0("ref.g", genes.values), 
		paste0("ref.sg", selg.values), 
		paste0("ref.d", signif(d.values, digits=2)),
		paste0("ref.s", signif(s.values, digits=2))))
	
	
################# Running simulations

list.sim <- mclapply(torun, function(ff) eval(ff)(), mc.cores=min(length(torun), ceiling(param$mc.cores/param$sim.reps)))

################# Figure

ww              <- c("fitness", "initenv", "lateenv","initmut","latemut", "stability")
default.values  <- c(m=default$mut.rate, sm=default$sim.mutsd, a=default$a, N=default$N, g=default$n, sg=default$nsel, d=default$density, s=default$s)

pdf("figS7.pdf", width=1.2*param$maxfigwidth/param$figscale, height=param$maxfigwidth/param$figscale, pointsize=0.8*param$pointsize)
	layout(matrix(1:(length(ww)*length(captions)), ncol=length(captions)))
	par(mar=c(0.5, 0.5, 0.5, 0.1), oma=c(5, 4, 0, 0), cex=1)
	for (pp in names(captions)) {
		for (what in ww) {
			ls <- list.sim[grep(names(list.sim), pattern=paste0("ref\\.",pp,"\\d"))]
			xval <- sapply(strsplit(names(ls), split=paste0("\\.",pp)), function(sp) as.numeric(sp[2]))
			yval <- sapply(ls, function(x) mean(x$mean[[as.character(param$G)]][[what]]))
			yvar <- sapply(ls, function(x) mean(x$var[[as.character(param$G)]][[what]]))
			
			plot(NULL, log=if(pp %in% c("g","sg","a","d")) "" else "x", xaxt="n", yaxt="n", xlab="", ylab="", xlim=range(xval), ylim=ylims[[what]])
			arrows(x0=xval, y0=yval-sqrt(yvar), y1=yval+sqrt(yvar), code=3, length=0, angle=90, col="darkgray")
			points(xval, yval, type="p")
			
			abline(v=default.values[pp], col="gray", lty=3, lwd=3)
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
