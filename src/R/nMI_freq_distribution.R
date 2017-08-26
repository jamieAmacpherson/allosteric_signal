#! /usr/bin/R

# Script calculates the frequency distribution and semi-log frequency distribution of a MI matrix

dif = as.matrix(read.table("diff_mat.dat"))
apo = as.matrix(read.table("apo_nMI.mat"))
fbp = as.matrix(read.table("fbp_nMI.mat"))


MIdistribution = function(nmi_matrix, difference){
	mathist = hist(nmi_matrix, plot=FALSE, breaks=2000)
	pfdr = mean(nmi_matrix) + (3*sd(nmi_matrix))
	nfdr = mean(nmi_matrix) - (3*sd(nmi_matrix))
	plot(mathist$breaks[-1],
		mathist$counts,
		type='h',
	#	log='y',
		xlab = "nMI",
		ylab = expression(italic(f(nMI))),
		cex.lab = 2,
		cex.axis = 2,
		panel.first = grid())
	
	# line at mu + 3sd 	
	abline(v=pfdr,
		lty=2,
		col="red",
		lwd=2)
	
	# line at mu - 3sd
	abline(v=nfdr,
		lty=2,
		col="red",
		lwd=2)

	# text
	text(pfdr + 0.2, 100000, paste("nMI > ", round(pfdr, 3)), cex=2)
	if(difference == TRUE){
		text(nfdr - 0.2, 100000, paste("nMI < ", round(nfdr, 3)), cex=2)
	}
}

diflevelplot = function(nmi_matrix, ncutoff, pcutoff){
	
	library(lattice)

	cols = colorRampPalette(c("white", "forestgreen"),
			space = "rgb")

	cutoffs = c(min(nmi_matrix),
		pcutoff,
		max(nmi_matrix))

	x.scale <- list(at=seq(1,ncol(nmi_matrix), 20))
	y.scale <- list(at=seq(1,nrow(nmi_matrix), 100))

	levelplot(dif,
		xlab="Fragment",
		ylab="Fragment",
		col.regions = cols(2),
		cuts=2,
		at=cutoffs,
		cex.lab=2,
		scales=list(x=x.scale, y=y.scale))	
}


#pdf("nMI_freq_dist.pdf")
#par(mar=c(5,5,1,1))
#par(mfrow=c(3,1))
#MIdistribution(apo, FALSE)
#MIdistribution(fbp, FALSE)
#MIdistribution(dif, TRUE)

#dev.off()


pdf("nMI_difplot.pdf")
diflevelplot(dif,-0.02, 0.03)

dev.off()

