# 'source' this file before running any figure script

# Computing options

library(parallel)
default.mc.cores <- min(detectCores()-1, 128)

# General parameters

default.mut.correlated <- TRUE
default.a              <- 0.2
default.n              <- 6
default.dev.steps      <- 16
default.dev.measure    <- 4
default.epsilon.zero   <- 0.01

# Random networks

default.density        <- 1
#~ default.rand.mean      <- 0.07
#~ default.rand.sd        <- 0.18
default.rand.mean      <- 0
default.rand.sd        <- 1

# Robustness index parameters

default.rob.reps      <- 100
default.initenv.sd    <- 0.1
default.lateenv.sd    <- 0.1
default.initmut.sd    <- 0.1
default.latemut.sd    <- 0.1
default.log.robustness<- TRUE

# Simulation parameters

default.sim.reps      <- 20
default.G             <- 100
default.summary.every <- 10
default.N             <- 1000
default.initsd        <- 0.001
default.mut.rate      <- 0.01
default.sim.mutsd     <- 0.1
default.s             <- 10
default.nsel          <- 3
