#!/usr/bin/env Rscript

#########################################################
# simul.R
#
# Simulates the Wagner model
# 
# Copyright Arnaud Le Rouzic / CNRS 2019
# <lerouzic@egce.cnrs-gif.fr>
#
# Released under the WTFPL version 2.0
# * No warranty *
###########################################################

library(parallel) # loads the package even if the user does not want to run it parallel

# get the model.M2 function (the genotype-phenotype routine)
library(funr)
curr_dir <- funr::get_script_path()
if (length(curr_dir)==0) curr_dir <- "."
suppressPackageStartupMessages(source(paste(curr_dir, "netw.R", sep="/")))


######## Helper functions
rtruncnorm <- function(N, mean = 0, sd = 1, a = -Inf, b = Inf) {
  if (a > b) stop('Error: Truncation range is empty');
  U <- runif(N, pnorm(a, mean, sd), pnorm(b, mean, sd));
  qnorm(U, mean, sd); }
  

######## Parameter handling
simulation.param <- list(
	n=5, # Size of the gene network
	a=0.2, # Constitutive expression
	r=0.5, # Recombination rate
	init.sd = 0.01,  # initial distribution of alleles
	mut.correlated = TRUE, # Whether the new allele is centered on the old one
	dev.steps=20,
	G=100, # Number of generations
	summary.every=1, # Number of generations between summary points
	debug = FALSE,
	outfile="out.txt",
	mc.cores = 1
)

generation.param <- list(
	N = 1000,
	mut.rate = 0.01, # per individual
	mut.sd = 0.1,
	som.rate = 0.,
	som.sd = 0.,
	initenv.sd = 0,
	lateenv.sd = 0.)
	

fitness.param <- list(
	theta = 0.5, 
	s = 10)

test.param <- list(
	test.rep = 1000,
	test.initenv.sd = 1,
	test.lateenv.sd = 0.1, 
	test.initmut.sd = 0.1,
	test.latemut.sd = 0.1,
	test.indiv = FALSE)

printhelp <- function() {
	fullist <- c(simulation.param, generation.param, fitness.param, test.param)
	cat("Printing help\n")
	cat("\tParameters and default values\n")
	for (i in seq_along(fullist)) {
		cat(paste0("\t\t-", names(fullist)[i], " = "))
		for (j in seq_along(fullist[[i]])) 
			cat(fullist[[i]][j], " ")
		cat("\n")
	}
	cat("\n")
}

parsecommandline <- function() {
	argss <- commandArgs(trailingOnly=TRUE)
	if (length(argss) == 0) return(list())
	whichpar <- grep(argss, pattern="^-")
	pars <- sapply(strsplit(argss[whichpar], split='-'), function(x) x[2])
	if (any (pars == "h")) {
		printhelp()
		quit()
	}
	nextpar <- c(whichpar, length(argss)+1)
	whichval <- lapply(whichpar, function(i) { upto <- nextpar[which(nextpar > i)[1]]-1; if (upto < i+1) numeric(0) else seq(i+1, upto)} )
	undef <- which(sapply(whichval, length) < 1)
	if(length(undef) > 0) stop("Parameter ", paste(pars[undef], collapse=", "), ": no value provided.")
	fulllist <- c(simulation.param, generation.param, fitness.param, test.param)
	unknown <- pars[!pars %in% names(fulllist)]
	if (length(unknown) > 0) {
		whichval[pars %in% unknown] <- NULL
		pars[pars %in% unknown] <- NULL
		warning("Parameter ", paste0(unknown, collapse=", "), " unknown and will be ignored.")
	}
	ans <- lapply(whichval, function(x) argss[x])
	names(ans) <- pars
	for (i in seq_along(ans)) mode(ans[[i]]) <- mode(fulllist[[names(ans)[i]]])
	ans
}



######### Output management 

summarypop <- function(pop, param.sim, param.test) {
	if (!is.list(pop[[1]])) browser()
	indiv.genot <- lapply(pop, function(i) 0.5*i$genotype.father + 0.5*i$genotype.mother)
	ans <- list(
		genotype.mean = apply(simplify2array(indiv.genot), 1:2, mean),
		genotype.var  = apply(simplify2array(indiv.genot), 1:2, var),
		phenotype.mean = colMeans(do.call(rbind, lapply(pop, function(x) x$phenotype[[1]]))),
		phenotype.var  = var(do.call(rbind, lapply(pop, function(x) x$phenotype[[1]]))),
		fitness.mean   = mean(sapply(pop, function(x) x$fitness)), 
		fitness.var    = var(sapply(pop, function(x) x$fitness))
	)
	if (param.test$test.indiv) {
		ans$robustness.initenv <- rowMeans(do.call(cbind, mclapply(indiv.genot, 
			function(W) robustness.initenv(W=W, param.sim=param.sim, env.sd=param.test$test.initenv.sd, rep=param.test$test.rep), mc.cores=param.sim$mc.cores)))
		ans$robustness.initenv <- rowMeans(do.call(cbind, lapply(indiv.genot, 
			function(W) robustness.lateenv(W=W, param.sim=param.sim, env.sd=param.test$test.lateenv.sd, rep=param.test$test.rep), mc.cores=param.sim$mc.cores)))
		ans$robustness.initenv <- rowMeans(do.call(cbind, lapply(indiv.genot, 
			function(W) robustness.initmut(W=W, param.sim=param.sim, mut.sd=param.test$test.initmut.sd, rep=param.test$test.rep), mc.cores=param.sim$mc.cores)))
		ans$robustness.initenv <- rowMeans(do.call(cbind, lapply(indiv.genot, 
			function(W) robustness.latemut(W=W, param.sim=param.sim, mut.sd=param.test$test.latemut.sd, rep=param.test$test.rep), mc.cores=param.sim$mc.cores)))
		ans$robustness.initenv <- rowMeans(do.call(cbind, lapply(indiv.genot, 
			function(W) robustness.stability(W=W, param.sim=param.sim), mc.cores=param.sim$mc.cores)))
	} else {
		Wmean <- matrix(ans$genotype.mean, ncol=sqrt(length(ans$genotype.mean)))
		ans$robustness.initenv <- robustness.initenv(W=Wmean, param.sim=param.sim, env.sd=param.test$test.initenv.sd, rep=param.test$test.rep)
		ans$robustness.lateenv <- robustness.lateenv(W=Wmean, param.sim=param.sim, env.sd=param.test$test.lateenv.sd, rep=param.test$test.rep)
		ans$robustness.initmut <- robustness.initmut(W=Wmean, param.sim=param.sim, mut.sd=param.test$test.initmut.sd, rep=param.test$test.rep)
		ans$robustness.latemut <- robustness.latemut(W=Wmean, param.sim=param.sim, mut.sd=param.test$test.latemut.sd, rep=param.test$test.rep)
		ans$robustness.stability <- robustness.stability(W=Wmean, param.sim=param.sim)
	}
	ans
}

summaryprint <- function(summarylist, file=NULL) {
	final.df <- do.call(rbind, lapply(summarylist, function(ss) {
		unlist(ss)
	}))
	rownames(final.df) <- names(summarylist)
	if (!is.null(file))
		write.table(final.df, file=file, sep="\t", quote=FALSE)
	invisible(final.df)
}


######## Checking the consistency of the state of the program

check.indiv <- function(indiv) {
	stopifnot(!is.null(indiv), !is.list(indiv))
	stopifnot(all(c("genotype.father", "genotype.mother", "phenotype", "fitness") %in% names(indiv)))
	stopifnot(is.matrix(indiv$genotype.father), is.numeric(indiv$genotype.father), sum(is.na(indiv$genotype.father)) == 0, 
		nrow(indiv$genotype.father) > 0, ncol(indiv$genotype.father) == nrow(indiv$genotype.father))
	stopifnot(is.matrix(indiv$genotype.mother), is.numeric(indiv$genotype.mother), sum(is.na(indiv$genotype.mother)) == 0, 
		nrow(indiv$genotype.mother) > 0, ncol(indiv$genotype.mother) == nrow(indiv$genotype.mother))
	stopifnot(nrow(indiv$genotype.father) == nrow(indiv$genotype.mother))
	stopifnot(sum(is.na(indiv$phenotype))==0, sum(indiv$phenotype < 0) == 0, sum(indiv$phenotype > 1) == 0)
	stopifnot(length(indiv$phenotype) == nrow(indiv$genotype.father))
	stopifnot(!is.na(indiv$fitness))
}

check.population <- function(population) {
	invisible(sapply(population, check.indiv))
}


########## Genotype-phenotype relationships

GPmap <- function(W, param.sim, param.gen) {
	if(!is.matrix(W)) W <- matrix(W, ncol=sqrt(length(W)))
	stopifnot(ncol(W) == nrow(W))
	S0 <- rtruncnorm(nrow(W), mean=param.sim$a, sd=param.gen$initenv.sd, 0, 1) 
	
	phen1 <- model.M2(W, a=param.sim$a, S0=S0, steps=param.sim$dev.steps, measure=1)
	
	if (param.gen$som.rate > 0 && param.gen$som.sd > 0) {
		nbmut <- min(rpois(1, param.gen$som.rate), sum(W != 0))
		if (nbmut > 0) {
			which.mut <- sample(size=nbmut,  which(W != 0), replace=FALSE) # Bug if only one W != 0
			mm <- if(param.sim$mut.correlated) W[which.mut] else 0
			W[which.mut] <- rnorm(nbmut, mean=mm, sd=param.gen$som.sd)
		}
	} 
	S02 <- rtruncnorm(nrow(W), phen1$mean, sd=param.gen$initenv.sd, 0, 1) 
	phen2 <- model.M2(W, a=param.sim$a, S0=S02, steps=1, measure=1)
	return(list(phen1=phen1$mean, phen2=phen2$mean))
}

fitness <- function(phenotype, param.fit) {
	ans <- prod(exp(-param.fit$s*(phenotype$phen2-param.fit$theta)^2))
}

########## Calculation of robustness scores

robustness.initenv <- function(W, param.sim, env.sd, rep=1000) {
	ans <- replicate(rep,  
		model.M2(W, a=param.sim$a, S0=rtruncnorm(nrow(W), mean=param.sim$a, sd=env.sd, 0, 1) , steps=param.sim$dev.steps, measure=1)$mean)
	apply(ans, 1, var)
}

robustness.lateenv <- function(W, param.sim, env.sd, rep=1000) {
	ref <- model.M2(W, a=param.sim$a, S0=rep(param.sim$a, ncol(W)) , steps=param.sim$dev.steps, measure=1)$mean
	ans <- replicate(rep, 
		model.M2(W, a=param.sim$a, S0=rtruncnorm(nrow(W), mean=ref, sd=env.sd, 0, 1) , steps=1, measure=1)$mean)
	apply(ans, 1, var)
}

robustness.initmut <- function(W, param.sim, mut.sd, nbmut=1, rep=1000) {
	ans <- replicate(rep,  
		{
			myW <- W
			nbmut <- min(nbmut, sum(W != 0))
			which.mut <- sample(size=nbmut,  which(W != 0), replace=FALSE) # Bug if only one W != 0
			mm <- if(param.sim$mut.correlated) W[which.mut] else 0
			myW[which.mut] <- rnorm(nbmut, mean=mm, sd=mut.sd)
			model.M2(myW, a=param.sim$a, S0=rep(param.sim$a, ncol(W)), steps=param.sim$dev.steps, measure=1)$mean
		})
	apply(ans, 1, var)
}

robustness.latemut <- function(W, param.sim, mut.sd, nbmut=1, rep=1000) {
	ref <- model.M2(W, a=param.sim$a, S0=rep(param.sim$a, ncol(W)) , steps=param.sim$dev.steps, measure=1)$mean
	ans <- replicate(rep,  
		{
			myW <- W
			nbmut <- min(nbmut, sum(W != 0))
			which.mut <- sample(size=nbmut,  which(W != 0), replace=FALSE) # Bug if only one W != 0
			mm <- if(param.sim$mut.correlated) W[which.mut] else 0
			myW[which.mut] <- rnorm(nbmut, mean=mm, sd=mut.sd)
			model.M2(myW, a=param.sim$a, S0=ref, steps=1, measure=1)$mean
		})
	apply(ans, 1, var)
}

robustness.stability <- function(W, param.sim) {
	ref <- model.M2(W, a=param.sim$a, S0=rep(param.sim$a, ncol(W)) , steps=param.sim$dev.steps, measure=1)$mean
	onemore <- model.M2(W, a=param.sim$a, S0=ref, steps=1, measure=1)$mean
	(ref-onemore)^2
}


########## Simulation routines

initpop <- function(param.sim, param.gen, param.fit) {
	.initindiv <- function() { 
		genotype.father <- matrix(rnorm(param.sim$n^2, mean=0, sd=param.sim$init.sd), ncol=param.sim$n)
		genotype.mother <- matrix(rnorm(param.sim$n^2, mean=0, sd=param.sim$init.sd), ncol=param.sim$n)
		phenotype <- GPmap(0.5*genotype.father + 0.5*genotype.mother, param.sim, param.gen)
		fitness <- fitness(phenotype, param.fit)
		list(genotype.father=genotype.father, genotype.mother=genotype.mother, phenotype=phenotype, fitness=fitness)
	}
	
	mclapply(1:param.gen$N, function(i) .initindiv(), mc.cores=param.sim$mc.cores)
}

makegamete <- function(indiv, param.sim, param.gen) {
	rec <- rbinom(param.sim$n, 1, prob=c(0.5, rep(param.sim$r, param.sim$n-1)))
	is.father <- cumsum(rec)%%2==0 # a bit complex, but should do the job fast
	ans <- do.call(rbind, lapply(seq_along(is.father), function(i) if (is.father[i]) indiv$genotype.father[i,] else indiv$genotype.mother[i,]))
	nbmut <- min(rpois(1, param.gen$mut.rate), sum(ans != 0))
	if (nbmut > 0) {
		which.mut <- sample(size=nbmut, which(ans != 0), replace=FALSE)
		mm <- if(param.sim$mut.correlated) ans[which.mut] else 0
		ans[which.mut] <- rnorm(nbmut, mean=mm, sd=param.gen$mut.sd)
	}
	matrix(ans, ncol=sqrt(length(ans)))
}

generation <- function(population, param.sim, param.gen, param.fit) {
	.newindiv <- function(genotype.father, genotype.mother) {
		phenotype <- GPmap(0.5*genotype.father+ 0.5*genotype.mother, param.sim, param.gen)
		fitness <- fitness(phenotype, param.fit)
		list(genotype.father=genotype.father, genotype.mother=genotype.mother, phenotype=phenotype, fitness=fitness)
	}
	
	if (param.sim$debug) check.population(population)
	
	# hermaphrodite population
	fitnesses <- sapply(population, "[", "fitness")
	fathers <- sample(size=param.gen$N, seq_along(population), prob=fitnesses, replace=TRUE)
	mothers <- sample(size=param.gen$N, seq_along(population), prob=fitnesses, replace=TRUE)
	
	newpop <- mcmapply(fathers, mothers, FUN=function(father.id, mother.id) .newindiv(makegamete(population[[father.id]], param.sim, param.gen), makegamete(population[[mother.id]], param.sim, param.gen)), SIMPLIFY=FALSE, mc.cores=param.sim$mc.cores)
	if (param.sim$debug) check.population(newpop)
	newpop
}

simulation <- function(param.sim, param.gen, param.fit, param.test) {
	summary.ans <- list()
	pop <- initpop(param.sim, param.gen, param.fit)
	for (gg in 1:param.sim$G) {
		if (gg %% param.sim$summary.every == 0) summary.ans[[as.character(gg)]] <- summarypop(pop, param.sim, param.test)
		pop <- generation(pop, param.sim, param.gen, param.fit)
	}
	summary.ans[[as.character(param.sim$G)]] <- summarypop(pop, param.sim, param.test)
	summary.ans
}

################ START OF THE EXECUTABLE SCRIPT #####################

#~ Rprof("prof.txt")

mypar <- parsecommandline()
simulation.param[names(simulation.param) %in% names(mypar)] <- mypar[names(mypar) %in% names(simulation.param)]
generation.param[names(generation.param) %in% names(mypar)] <- mypar[names(mypar) %in% names(generation.param)]
fitness.param[names(fitness.param) %in% names(mypar)] <- mypar[names(mypar) %in% names(fitness.param)]
test.param[names(test.param) %in% names(mypar)] <- mypar[names(mypar) %in% names(test.param)]

mysim <- simulation(simulation.param, generation.param, fitness.param, test.param)
summaryprint(mysim, file=simulation.param$outfile)

#~ Rprof(NULL)
