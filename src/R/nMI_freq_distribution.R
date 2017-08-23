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
		log='y',
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

pdf("nMI_freq_logdist.pdf")
par(mar=c(5,5,1,1))
par(mfrow=c(3,1))
MIdistribution(apo, FALSE)
MIdistribution(fbp, FALSE)
MIdistribution(dif, TRUE)

dev.off()
