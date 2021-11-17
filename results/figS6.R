#!/usr/bin/env Rscript

# Illustrates the robustness differences between all study cases

source("./studycases.R")
source("./terminology.R")
source("./defaults.R")

source("../src/netw.R")
source("../src/tools.R")

################ Options

param <- default

cols <- c("black","red")

################ Graphical functions

illustrate.reference <- function(W, ...) {
	mod <- model.M2(W=W, a=param$a, steps=param$dev.steps, full=TRUE)
	plot(NULL, xlim=c(1, param$dev.steps+2), ylim=c(0,1), xaxt="n", yaxt="n", xlab="", ylab="")
	for (i in 1:nrow(W))
		lines(mod$full[i,], col=cols[i], ...)
}

illustrate.initmut <- function(W, rep=10, ...) {
	for (r in 1:rep) {
		myW <- W
		whichmut <- sample(seq_along(W), 1)
		myW[whichmut] <- rnorm(1, myW[whichmut], sd=param$initmut.sd)
		mod <- model.M2(W=myW, a=param$a, steps=param$dev.steps, full=TRUE)
		for (i in 1:nrow(W))
			lines(mod$full[i,], col=makeTransparent(cols[i]), ...)
	}
}

illustrate.initenv <- function(W, rep=10, ...) {
	for (r in 1:rep) {
		myS0 <- rnorm(nrow(W), param$a, sd=param$initenv.sd)
		myS0[myS0<0] <- 0
		myS0[myS0>1] <- 1
		mod <- model.M2(W=W, a=param$a, S0=myS0, steps=param$dev.steps, full=TRUE)
		for (i in 1:nrow(W))
			lines(mod$full[i,], col=makeTransparent(cols[i]), ...)
	}
}

illustrate.latemut <- function(W, rep=20, ...) {
	ref <- model.M2(W=W, a=param$a, steps=param$dev.steps)
	for (r in 1:rep) {
		myW <- W
		whichmut <- sample(seq_along(W), 1)
		myW[whichmut] <- rnorm(1, myW[whichmut], sd=param$latemut.sd)
		mod <- model.M2(W=myW, a=param$a, S0=ref$mean, steps=1, full=TRUE)
		for (i in 1:nrow(W))
			lines(x=(1:2)+param$dev.steps, mod$full[i,], col=makeTransparent(cols[i]), ...)
	}
}

illustrate.lateenv <- function(W, rep=20, ...) {
	ref <- model.M2(W=W, a=param$a, steps=param$dev.steps)
	for (r in 1:rep) {
		myS0 <- rnorm(nrow(W), ref$mean, sd=param$lateenv.sd)
		myS0[myS0<0] <- 0
		myS0[myS0>1] <- 1
		mod <- model.M2(W=W, a=param$a, S0=myS0, steps=1, full=TRUE)
		for (i in 1:nrow(W))
			lines(x=(1:2)+param$dev.steps, mod$full[i,], col=makeTransparent(cols[i]), ...)
	}
}


#################### Figure

pdf("figS6.pdf", width=param$maxfigwidth/param$figscale, height=0.9*param$maxfigwidth/param$figscale, pointsize=param$pointsize)
	layout(matrix(1:(4*nrow(stud)), ncol=4, byrow=TRUE))
	par(mar=c(0,0,0,0)+0.2, oma=c(4, 4, 4, 0), xpd=NA, cex=1)
	
	for (rstud in 1:nrow(stud)) {
		# target is defined in studycases.R
		W <- targetW(cbind(stud[rstud,], rep(NA, network.size)), target=target, a=param$a) 
		
		illustrate.reference(W, lwd=2)
		illustrate.initmut(W)
		axis(2)
		mtext("Expression", 2, outer=FALSE, line=3)
		if (rstud == 1) 
			mtext(TERM.GENCAN.LONG, 3, line=2, cex=1.5)
		if (rstud == nrow(stud)) {
			axis(1)
			mtext("Time steps", 1, outer=FALSE, line=3)
		}
		subpanel(LETTERS[rstud], cex=1.8, line=-1.5)
		
		illustrate.reference(W, lwd=2)
		illustrate.latemut(W)
		if (rstud == 1) 
			mtext(TERM.SOM.LONG, 3, line=2, cex=1.5)	
		if (rstud == nrow(stud)) {
			axis(1)
			mtext("Time steps", 1, outer=FALSE, line=3)
		}
			
		illustrate.reference(W, lwd=2)
		illustrate.initenv(W)
		if (rstud == 1) 
			mtext(TERM.ENVCAN.LONG, 3, line=2, cex=1.5)
		if (rstud == nrow(stud)) {
			axis(1)
			mtext("Time steps", 1, outer=FALSE, line=3)
		}
			
		illustrate.reference(W, lwd=2)
		illustrate.lateenv(W)
		if (rstud == 1) 
			mtext(TERM.HOMEO.LONG, 3, line=2, cex=1.5)
		if (rstud == nrow(stud)) {
			axis(1)
			mtext("Time steps", 1, outer=FALSE, line=3)
		}		
	}
dev.off()
