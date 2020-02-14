source("commonpure.R")

n.genes <- 6
sel.genes <- 3
s <- c(rep(10, sel.genes), rep(0, n.genes-sel.genes))
W0 <- matrix(rnorm(n.genes^2, sd=0.000001), ncol=n.genes)
reps <- 3
test.rep <- 10
grad.effect <- 0.1
N <- 500
G <- 1000
force.run <- TRUE


# Null model: no direct selection on robustness
null.sim <- pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=10, grad.rob=c(0,0,0,0,0)), reps=reps, series.name="pure-null", force.run=force.run)

ie.m <- pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=10, grad.rob=c(-grad.effect,0,0,0,0)), reps=reps, series.name="pure-ie-m", force.run=force.run)
le.m <- pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=10, grad.rob=c(0,-grad.effect,0,0,0)), reps=reps, series.name="pure-le-m", force.run=force.run)
im.m <- pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=10, grad.rob=c(0,0,-grad.effect,0,0)), reps=reps, series.name="pure-im-m", force.run=force.run)
lm.m <- pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=10, grad.rob=c(0,0,0,-grad.effect,0)), reps=reps, series.name="pure-lm-m", force.run=force.run)
st.m <- pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=10, grad.rob=c(0,0,0,0,-grad.effect)), reps=reps, series.name="pure-st-m", force.run=force.run)

ie.p <- pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=10, grad.rob=c(grad.effect,0,0,0,0)), reps=reps, series.name="pure-ie-p", force.run=force.run)
le.p <- pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=10, grad.rob=c(0,grad.effect,0,0,0)), reps=reps, series.name="pure-le-p", force.run=force.run)
im.p <- pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=10, grad.rob=c(0,0,grad.effect,0,0)), reps=reps, series.name="pure-im-p", force.run=force.run)
lm.p <- pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=10, grad.rob=c(0,0,0,grad.effect,0)), reps=reps, series.name="pure-lm-p", force.run=force.run)
st.p <- pure.run.reps(W0, list(s=s, G=G, N=N, rep=test.rep, summary.every=10, grad.rob=c(0,0,0,0,grad.effect)), reps=reps, series.name="pure-st-p", force.run=force.run)

plot(as.numeric(names(null.sim$mean)), sapply(null.sim$mean, "[", "fitness"), type="b")
