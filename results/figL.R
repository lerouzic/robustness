#!/usr/bin/env Rscript

source("./commonpure.R")
source("./terminology.R")

library(parallel)
mc.cores <- min(detectCores()-1, 128)

phen <- c(
	fitness="Fitness",
    initenv=substitute(x~(y), list(x=TERM.ENVCAN.LONG, y=ABBRV.ENVCAN[[1]])),
    lateenv=substitute(x~(y), list(x=TERM.HOMEO.LONG, y=ABBRV.HOMEO[[1]])),
    initmut=substitute(x~(y), list(x=TERM.GENCAN.LONG, y=ABBRV.GENCAN[[1]])),
    latemut=substitute(x~(y), list(x=TERM.SOM.LONG, y=ABBRV.SOM[[1]])),
    stability=substitute(x~(y), list(x=TERM.STAB.LONG, y=ABBRV.STAB[[1]])))
    

n.genes <- 6
sel.genes <- 3
s <- c(rep(10, sel.genes), rep(0, n.genes-sel.genes))
W0 <- matrix(rnorm(n.genes^2, sd=0.000001), ncol=n.genes)
reps <- 20
test.rep <- 10
N <- 1000
G <- 5000
force.run <- FALSE
mut.rate <- 0.001
nb.values <- 11

mut.values <- 10^seq(-6,-1, length.out=nb.values)
N.values   <- round(10^seq(1, 4, length.out=nb.values))
genes.values <- round(seq(3, 20, length.out=nb.values))
selg.values  <- 1:6
s.values <- 10^seq(-2, 4, length.out=nb.values)

torun.mut <- lapply(mut.values, function(mm) 
	substitute(function() pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=G, mut.rate=mm), 
		reps=reps, series.name=paste0("real-ref-m", signif(mm, digits=2)), force.run=force.run), list(mm=mm)))
torun.N <- lapply(N.values, function(nn) 
	substitute(function() pure.run.reps(W0, list(s=s, G=G, N=nn, rep=test.rep, summary.every=G, mut.rate=mut.rate), 
		reps=reps, series.name=paste0("real-ref-N", nn), force.run=force.run), list(nn=nn)))
torun.genes <- lapply(genes.values, function(gg) 
	substitute(function() pure.run.reps(W0=matrix(rnorm(gg^2, sd=0.000001), ncol=gg), list(s=c(rep(10, sel.genes), rep(0 ,gg-sel.genes)), G=G, N=N, rep=test.rep, summary.every=G, mut.rate=mut.rate), 
		reps=reps, series.name=paste0("real-ref-g", gg), force.run=force.run), list(gg=gg)))
torun.selg <- lapply(selg.values, function(sg)
	substitute(function() pure.run.reps(W0, list(s=c(rep(10,sg),rep(0,n.genes-sg)), G=G, N=N, rep=test.rep, summary.every=G, mut.rate=mut.rate), 
		reps=reps, series.name=paste0("real-ref-sg", sg), force.run=force.run), list(sg=sg)))
torun.sel <- lapply(s.values, function(ss)
	substitute(function() pure.run.reps(W0, list(s=c(rep(ss,sel.genes),rep(0,n.genes-sel.genes)), G=G, N=N, rep=test.rep, summary.every=G, mut.rate=mut.rate), 
		reps=reps, series.name=paste0("real-ref-s", signif(ss, digits=2)), force.run=force.run), list(ss=ss)))

torun <- setNames(
	c(torun.mut, torun.N, torun.genes, torun.selg, torun.sel),
	c(paste0("ref.m", signif(mut.values, digits=2)), paste0("ref.N", N.values), paste0("ref.g", genes.values), paste0("ref.sg", selg.values), paste0("ref.s", signif(s.values, digits=2))))
list.sim <- mclapply(torun, function(ff) eval(ff)(), mc.cores=min(length(torun), ceiling(mc.cores/reps)))

layout(matrix(1:25, ncol=5))
for (par in c("m","N","g","sg","s")) {
	for (what in c("initenv", "lateenv","initmut","latemut", "stability")) {
		ls <- list.sim[grep(names(list.sim), pattern=paste0("ref\\.",par,"\\d"))]
		val <- sapply(strsplit(names(ls), split=paste0("\\.",par)), function(sp) as.numeric(sp[2]))
		plot(val, sapply(ls, function(x) mean(x$mean[[as.character(G)]][[what]])), xlab=par, ylab=what, type="o", log=if(par %in% c("g","sg")) "" else "x")
	}
}
