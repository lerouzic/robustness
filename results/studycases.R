network.size <- 2
target <- c(0.3, 0.6)

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

# Study cases
stud <- rbind(
    A=c(0.7,0.2),
    B=c(-0.3,0.3),
    C=c(-0.4,0.8),
    D=c(-1,-0.8),
    E=c(1.5,3.5))
