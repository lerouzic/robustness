#!/usr/bin/env Rscript

# Illustrates the robustness differences between all study cases

source("../src/netw.R")
source("./studycases.R")
source("./terminology.R")
source("./defaults.R")

################ Options
cols <- 1:2

a          <- default.a
initmut.sd <- default.initmut.sd
latemut.sd <- default.latemut.sd
initenv.sd <- default.initenv.sd
lateenv.sd <- default.lateenv.sd
dev.steps <-  default.dev.steps

############## Functions
makeTransparent<-function(someColor, alpha=70)
{ # from https://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}


illustrate.reference <- function(W, ...) {
	mod <- model.M2(W=W, a=a, steps=dev.steps, full=TRUE)
	plot(NULL, xlim=c(1, dev.steps+2), ylim=c(0,1), xaxt="n", yaxt="n", xlab="", ylab="")
	for (i in 1:nrow(W))
		lines(mod$full[i,], col=cols[i], ...)
}

illustrate.initmut <- function(W, rep=10, ...) {
	for (r in 1:rep) {
		myW <- W
		whichmut <- sample(seq_along(W), 1)
		myW[whichmut] <- rnorm(1, myW[whichmut], sd=initmut.sd)
		mod <- model.M2(W=myW, a=a, steps=dev.steps, full=TRUE)
		for (i in 1:nrow(W))
			lines(mod$full[i,], col=makeTransparent(cols[i]), ...)		
	}
}

illustrate.initenv <- function(W, rep=10, ...) {
	for (r in 1:rep) {
		myS0 <- rnorm(nrow(W), a, sd=initenv.sd)
		myS0[myS0<0] <- 0
		myS0[myS0>1] <- 1
		mod <- model.M2(W=W, a=a, S0=myS0, steps=dev.steps, full=TRUE)
		for (i in 1:nrow(W))
			lines(mod$full[i,], col=makeTransparent(cols[i]), ...)		
	}
}

illustrate.latemut <- function(W, rep=20, ...) {
	ref <- model.M2(W=W, a=a, steps=dev.steps)
	for (r in 1:rep) {
		myW <- W
		whichmut <- sample(seq_along(W), 1)
		myW[whichmut] <- rnorm(1, myW[whichmut], sd=latemut.sd)
		mod <- model.M2(W=myW, a=a, S0=ref$mean, steps=1, full=TRUE)
		for (i in 1:nrow(W))
			lines(x=(1:2)+dev.steps, mod$full[i,], col=makeTransparent(cols[i]), ...)		
	}
}

illustrate.lateenv <- function(W, rep=20, ...) {
	ref <- model.M2(W=W, a=a, steps=dev.steps)
	for (r in 1:rep) {
		myS0 <- rnorm(nrow(W), ref$mean, sd=lateenv.sd)
		myS0[myS0<0] <- 0
		myS0[myS0>1] <- 1
		mod <- model.M2(W=W, a=a, S0=myS0, steps=1, full=TRUE)
		for (i in 1:nrow(W))
			lines(x=(1:2)+dev.steps, mod$full[i,], col=makeTransparent(cols[i]), ...)		
	}
}


#################### Figure

pdf("figS5.pdf", width=12, height=12)
	layout(matrix(1:(4*nrow(stud)), ncol=4, byrow=TRUE))
	par(mar=c(0,0,0,0)+0.2, oma=c(4, 4, 4, 0), xpd=NA)
	
	for (rstud in 1:nrow(stud)) {
		W <- targetW(cbind(stud[rstud,], rep(NA, network.size)), target=target, a=a)
		
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
		text(4, 0.95, LETTERS[rstud], cex=2)
		
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
