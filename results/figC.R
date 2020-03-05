#!/usr/bin/env Rscript

library(parallel)
mc.cores <- min(64, detectCores()-1)

source("./terminology.R")
source("./studycases.R")
source("../src/robindex.R")

phen <- c(
    #mean=TERM.EXPRESSION,
    initenv=TERM.ENVCAN.LONG,
    lateenv=TERM.HOMEO.LONG,
    initmut=TERM.GENCAN.LONG,
    latemut=TERM.SOM.LONG,
    stability=TERM.STAB.LONG)

a <- 0.2
dev.steps <- 20

rob.reps <- 1000
rob.initenv.sd <- 0.1
rob.lateenv.sd  <- 0.1
rob.mut.sd     <- 0.1

targetW <- function(W, target, a) {
    lambda <- (1-a)/a
    mu <- 1/(a*(1-a))
    frev <- function(x) -log((1-x)/lambda/x) / mu
    stopifnot(nrow(W) == ncol(W), nrow(W) == length(target))
    stopifnot(all(apply(W, 1, function(x) sum(is.na(x))) == 1))
    stopifnot(all(target > 0), all(target < 1))
    ans <- W
    for (i in 1:nrow(W)) {
        miniT <- frev(target[i])
        whichmiss <- which(is.na(W[i,]))
        minians <- (miniT - sum(W[i, -whichmiss]*target[-whichmiss]))/target[whichmiss]
        ans[i, whichmiss] <- minians
    }
    ans
}

difftarget <- function(res, target) {
    df <- t(res[,grep("mean.", colnames(res))])-target
    apply(abs(df), 2, sum)
}

plotres <- function(res, crit="mean", stud=NULL, mask=NULL, contour=FALSE, mx = 0.2) {
    z <- res[,grep(colnames(res), pattern=crit)[1]] # take only the first gene
    z[z>mx] <- mx
    if (!is.null(mask)) z[mask] <- NA
    
    if (crit %in% names(phen)) main <- phen[crit] else main <- crit
    
    image(x=ww1, y=ww2, z=matrix(z, nrow=sqrt(nrow(res))), main=main, xlab=expression(W[11]), ylab=expression(W[21]), col=heat.colors(128)[1:100])
    if (contour) contour(x=ww1, y=ww2, z=matrix(z, nrow=sqrt(nrow(res))), add=TRUE) 
    if (!is.null(stud))
        invisible(sapply(rownames(stud), function(rn) text(x=stud[rn,1], y=stud[rn,2], rn, col="blue")))
}

illustrcase <- function(w11, w21, left=FALSE, right=FALSE, label="") {
    w <- targetW(cbind(c(w11, w21), rep(NA, sz)), target=target)
    xlim <- c(1, 25)
    plot.test.homeo(w, a=0.2, xlim=xlim, replicates=40, exagg=4)
    mtext(label, 3)
    if (left) { axis(2); mtext(TERM.EXPRESSION, 2, line=3, cex=0.8) }
    if (right) mtext(TERM.HOMEO.SHORT, 4, line=1.5)
    
    plot.test.envcan(w, a=0.2, xlim=xlim)
    if (left) { axis(2); mtext(TERM.EXPRESSION, 2, line=3, cex=0.8) }
    if (right) mtext(TERM.ENVCAN.SHORT, 4, line=1.5)

    plot.test.gencan(w, a=0.2, xlim=xlim)
    if (left) { axis(2); mtext(TERM.EXPRESSION, 2, line=3, cex=0.8) }
    if (right) mtext(TERM.GENCAN.SHORT, 4, line=1.5)    
    
    plot.test.somcan(w, a=0.2, xlim=xlim)
    if (left) { axis(2); mtext(TERM.EXPRESSION, 2, line=3, cex=0.8) }
    if (right) mtext(TERM.SOM.SHORT, 4, line=1.5) 
    
    axis(1)
    mtext(TERM.TIMESTEPS, 1, line=3, cex=0.8)
}

network.size <- 2
target <- c(0.3, 0.6)
grid.size <- 21

ww1 <- seq(-1.5, 2, length.out=grid.size)
ww2 <- seq(-1, 4, length.out=grid.size)

res <- expand.grid(ww1, ww2)

# Pattern 
#   a   NA
#   b   NA
# (should not be important)

res <- cbind(res, t(sapply(1:nrow(res), function(i) {
        w <- targetW(W=cbind(unlist(res[i,]), rep(NA, network.size)), target=target, a=a)
        return( c(w)[(1+network.size*(network.size-1)):(network.size*network.size)])
    })))
colnames(res) <- paste("W", outer(1:network.size, 1:network.size, paste, sep="."), sep=".")
res <- cbind(res, do.call(rbind, mclapply(1:nrow(res), function(i) {
        w <- matrix(unlist(res[i,]), nrow=network.size)
        c(mean=model.M2(w, a, steps=dev.steps)$mean,
          initenv=robindex.initenv(w, a, dev.steps, rob.initenv.sd, rep=rob.reps),
          lateenv=robindex.lateenv(w, a, dev.steps, rob.initenv.sd, rep=rob.reps),
          initmut=robindex.initmut(w, a, dev.steps, rob.mut.sd, rep=rob.reps),
          latemut=robindex.latemut(w, a, dev.steps, rob.mut.sd, rep=rob.reps),
          stability=robindex.stability(w, a, dev.steps))
    }, mc.cores=mc.cores)))


pdf("figC.pdf", width=12, height=8)
layout(rbind(1:3, c(4:5, 0)))
mm <- difftarget(res, target) > 0.15
for (ppp in names(phen))
	plotres(res, ppp, stud, mask=mm, contour=TRUE)
dev.off()


#~ layout(matrix(1:(4*nrow(stud)), nrow=4))
#~ par(mar=c(0.5,0.5,0.1,0.1), oma=c(4, 4, 2, 3))
#~ invisible(sapply(rownames(stud), function(rn) {
#~         illustrcase(stud[rn,1], stud[rn,2], left=(rn==rownames(stud)[1]), right=(rn==rownames(stud)[nrow(stud)]), label=rn)
#~     }))
    
#~ nicetab <- t(apply(stud, 1, function(x) targetW(cbind(x, rep(NA, length(x))), target)))
#~ colnames(nicetab) <- paste0("$W_{", c(11,21,12,22),"}$")
#~ library(xtable)
#~ print(xtable(nicetab), sanitize.colnames.function=identity)

