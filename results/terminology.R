TERM.EXPRESSION <- "Expression"
TERM.TIMESTEPS <- "Time steps"

TERM.STAB.SHORT <- "Stab."
TERM.STAB.LONG  <- "Stability"
TERM.HOMEO.SHORT <- "Late Env."
TERM.HOMEO.LONG  <- "Late Environmental"
TERM.ENVCAN.SHORT <- "Early Env."
TERM.ENVCAN.LONG  <- "Early Environmental"
TERM.GENCAN.SHORT <- "Early Gen."
TERM.GENCAN.LONG <- "Early Genetic"
TERM.SOM.SHORT   <- "Late Gen."
TERM.SOM.LONG    <- "Late Genetic"

ABBRV.STAB <- expression(rho[S])
ABBRV.HOMEO <- expression(rho[e])
ABBRV.ENVCAN <- expression(rho[E])
ABBRV.GENCAN <- expression(rho[M])
ABBRV.SOM <- expression(rho[m])

SDLETTER.HOMEO  <- "e"
SDLETTER.ENVCAN <- "E"
SDLETTER.GENCAN <- "M"
SDLETTER.SOM    <- "m"

COL.STAB <- "blue"
COL.HOMEO <- "lightgreen"
COL.ENVCAN <- "darkgreen"
COL.GENCAN <- "darkred"
COL.SOM <- "red"

phen.expression <- c(
    initenv=substitute(x~(y), list(x=TERM.ENVCAN.LONG, y=ABBRV.ENVCAN[[1]])),
    lateenv=substitute(x~(y), list(x=TERM.HOMEO.LONG, y=ABBRV.HOMEO[[1]])),
    initmut=substitute(x~(y), list(x=TERM.GENCAN.LONG, y=ABBRV.GENCAN[[1]])),
    latemut=substitute(x~(y), list(x=TERM.SOM.LONG, y=ABBRV.SOM[[1]])),
    stability=substitute(x~(y), list(x=TERM.STAB.LONG, y=ABBRV.STAB[[1]]))
)

default.labels <- c(
	initenv=ABBRV.ENVCAN, 
	lateenv=ABBRV.HOMEO, 
	initmut=ABBRV.GENCAN, 
	latemut=ABBRV.SOM, 
	stability=ABBRV.STAB)

default.cols   <- c(
	initenv=COL.ENVCAN, 
	lateenv=COL.HOMEO, 
	initmut=COL.GENCAN, 
	latemut=COL.SOM, 
	stability=COL.STAB)
