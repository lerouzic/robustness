source("./randnetwork.R")
source("./defaults.R")

cache.dir            <- "../cache"
evolved.file.pattern <- 'figG-null-\\d+.rds'
evolved.gen          <- NA    # NA: last generation of the simulations

epsilon.zero         <- default.epsilon.zero

evolved.files <- list.files(path=cache.dir, pattern=evolved.file.pattern, full.names=TRUE)

pdf("figN.pdf", width=15, height=4)

	# Kinetics
	gens <- c(100, 1000, 2000, 5000, 10000)
	xlim <- c(-0.2, 0.5)
	layout(t(seq_along(gens)))
	for (gg in gens) {
		Wlist <- Wlist.fromfiles(evolved.files, gen=as.character(gg))
		hh <- hist(unlist(Wlist), breaks=seq(-1, 1, by=0.01), plot=FALSE)
		hh$counts <- hh$counts[hh$breaks[-length(hh$breaks)] > xlim[1] & hh$breaks[-1] < xlim[2]]
		hh$counts[hh$counts == 0] <- NA
		hh$breaks <- hh$breaks[hh$breaks > xlim[1] & hh$breaks < xlim[2]]
		# gap.barplot(hh$counts/sum(hh$counts), gap=c(0.03, 0.3) ,ylim=c(0,0.3))
		bb <- barplot(hh$counts/sum(hh$counts, na.rm=TRUE), log="y", ylim=c(1e-4, 0.5), 
			ylab=if(gg==gens[1]) "Frequency (log scale)" else "", xlab=expression(W[ij]), xaxt="n", yaxt=if(gg==gens[1]) "s" else "n",
			col=ifelse(hh$breaks[-length(hh$breaks)] <= - epsilon.zero - 1e-6 | hh$breaks[-1] >= epsilon.zero + 1e-6, "gray", "red") )
		title(main=paste0("Generation ", gg, ", Density: ", round(mean(abs(unlist(Wlist)) <= epsilon.zero), digits=2)))
		rr   <- diff(range(bb))/diff(range(hh$breaks))
		tt   <- pretty(hh$breaks)
		atbb <- bb[1] + (tt-tt[1])*rr
		axis(1, at=atbb, labels=as.character(tt))
	}
	
dev.off()
