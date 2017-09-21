#! /usr/bin/R

# Script calculates the frequency distribution and semi-log frequency distribution of a MI matrix

dif = as.matrix(read.table("diff_mat.dat"))
apo = as.matrix(read.table("apo_nMI.mat"))
fbp = as.matrix(read.table("fbp_nMI.mat"))


MIdistribution = function(nmi_matrix, difference){
	mathist = hist(nmi_matrix, plot=FALSE, breaks=2000)
	pfdr = mean(nmi_matrix) + (2*sd(nmi_matrix))
	nfdr = mean(nmi_matrix) - (2*sd(nmi_matrix))
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
	text(pfdr + 0.3, 100000, paste("nMI > ", round(pfdr, 3)), cex=2)
	if(difference == TRUE){
		text(nfdr - 0.3, 100000, paste("nMI < ", round(nfdr, 3)), cex=2)
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

diflevelimage = function(nmi_matrix, ncutoff, pcutoff){
	
	# set all values inside of 2 sigmas of the mean NA
	adjmat = replace(nmi_matrix, nmi_matrix > ncutoff & nmi_matrix < pcutoff, NA)
	
	frags = seq(from=1, to=ncol(adjmat))	

	# plot the filtered matrix
	pdf("image_plt.pdf")
	image.plot(frags, frags, adjmat)
	dev.off()
	
}

circosplt = function(nmi_matrix, pcutoff, rangemin, rangemax){
	library(fields)
	library(circlize)	
	library(reshape2)

	# rename rows and columns
	if(ncol(nmi_matrix) == 515) {
		row.names(nmi_matrix) = seq(from=1, to=ncol(nmi_matrix))
		colnames(nmi_matrix) = seq(from=1, to=ncol(nmi_matrix))
	}
		else {
			row.names(nmi_matrix) = c(sprintf("A%s",seq(1:515)),
						sprintf("B%s",seq(516:1030)),
						sprintf("C%s",seq(1031:1545)),
						sprintf("D%s",seq(1546:2060)))  
		}
	
	# extract the upper triangle of the MI matrix
	nmi_matrix[lower.tri(nmi_matrix)] = 0	
	
	# extract desired range within the matrix
	nmi_matrix = nmi_matrix[rangemin:rangemax,rangemin:rangemax] 
	
	
	# set all values smaller than user-defined cutoff to zero
	print(pcutoff)
	posmat = replace(nmi_matrix, nmi_matrix < pcutoff, 0)
	
	# melt the matrix into a datafrane
	mltdat = melt(as.matrix(posmat))
	var2 = as.numeric(gsub("[^0-9]", "", mltdat$Var2))

	circdat = as.data.frame(cbind(mltdat$Var1,
				var2,
				as.numeric(round(mltdat$value, 3))))	

	# remove all correlations which have a zero mutual information
	row_sub = apply(circdat, 1, function(row) all(row !=0 ))
	circdatfilt = circdat[row_sub,]

	names(circdatfilt) = c("Var1", "Var2", "value")
	
	# order the dataframe
	circdat_ord = circdatfilt[order(circdatfilt$Var1, circdatfilt$Var2),]
	rownames(circdat_ord) = NULL

#	par(mar=c(5,5,2,2))

	chordDiagram(x = circdat_ord, annotationTrack = "grid", preAllocateTracks = 1,
		annotationTrackHeight = c(0.05, 0.1),
		transparency = 0.25,
		direction.type = c("diffHeight"), diffHeight  = -0.04,
		directional=1,
	)
	
	circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  		xlim = get.cell.meta.data("xlim")
  		ylim = get.cell.meta.data("ylim")
  		sector.name = get.cell.meta.data("sector.index")
  		circos.text(mean(xlim), ylim[1] + .1,
			sector.name,
			facing = "clockwise",
			niceFacing = TRUE,
			adj = c(0, 0.5),
			col = "lightgray")
		},
	bg.border = NA)

	# plot total mutual information content per residue
#	totcont = apply(posmat, 2, sum)
#	plot(totcont, type='h',
#		cex.axis = 2,
#		cex.lab = 2,
#		cex = 2,
#		xlab = "Fragment",
#		ylab = "Mutual information content")
}

pdf("chainA_circ.pdf")
circosplt(dif, 0.998, 1, 515)
dev.off()

pdf("chainB_circ.pdf")
circosplt(dif, 0.998, 516, 1030)
dev.off()

pdf("chainC_circ.pdf")
circosplt(dif, 0.998, 1031, 1545)
dev.off()

pdf("chainD_circ.pdf")
circosplt(dif, 0.998, 1546, 2060)
dev.off()


protav = (dif[1:515, 1:515] + dif[516:1030, 516:1030] + dif[1031:1545, 1031:1545] + dif[1546:2060, 1546:2060]) / 4

pdf("av_prot_circ.pdf") 
circosplt(protav, 0.998, 1, 515)
dev.off()

splitmatrix = function(difference_matrix, ncutoff, pcutoff){

	positive_diffMI = replace(difference_matrix, difference_matrix < pcutoff, 0)
	negative_diffMI = -1 * (replace(difference_matrix, difference_matrix > ncutoff, 0))

	write.table(positive_diffMI,
			file="positive_diffMI.out",
			sep=" ",
			col.names=F,
			row.names=F)
	
	write.table(negative_diffMI,
			file="negative_diffMI.out",
			sep=" ",
			col.names=F,
			row.names=F)
	

}


pdf("nMI_freq_dist.pdf")
par(mar=c(5,5,1,1))
par(mfrow=c(3,1))
MIdistribution(apo, FALSE)
MIdistribution(fbp, FALSE)
MIdistribution(dif, TRUE)

dev.off()


diflevelimage(dif,-0.046, 0.041)

splitmatrix(dif, -0.046, 0.041)

interprotcircosplt = function(nmi_matrix, pcutoff, xrangemin, xrangemax, yrangemin, yrangemax, outprefix){
	library(fields)
	library(circlize)	
	library(reshape2)

	# rename rows and columns
	row.names(nmi_matrix) = c(sprintf("A%s",seq(1:515)),
				sprintf("B%s",seq(516:1030)),
				sprintf("C%s",seq(1031:1545)),
				sprintf("D%s",seq(1546:2060)))  
	
	colnames(nmi_matrix) = c(sprintf("A%s",seq(1:515)),
				sprintf("B%s",seq(516:1030)),
				sprintf("C%s",seq(1031:1545)),
				sprintf("D%s",seq(1546:2060))) 
	
	# extract the upper triangle of the MI matrix
	nmi_matrix[lower.tri(nmi_matrix)] = 0	
	
	# extract desired range within the matrix
	nmi_matrix = nmi_matrix[xrangemin:xrangemax, yrangemin:yrangemax] 
	
	
	# set all values smaller than user-defined cutoff to zero
	print(pcutoff)
	posmat = replace(nmi_matrix, nmi_matrix < pcutoff, 0)
	
	# melt the matrix into a datafrane
	mltdat = melt(as.matrix(posmat))

	# remove all correlations which have a zero mutual information
	row_sub = apply(mltdat, 1, function(row) all(row !=0 ))
	circdatfilt = mltdat[row_sub,]

	
	# order the dataframe
	circdat_ord = mltdat[order(mltdat$Var1, mltdat$Var2),]
	rownames(circdat_ord) = NULL
	print(max(circdat_ord$value))

	par(mar=c(5,5,2,2))

	chordDiagram(x = circdat_ord, annotationTrack = "grid", preAllocateTracks = 1,
		annotationTrackHeight = c(0.05, 0.1),
		transparency = 0.25,
		direction.type = c("diffHeight"), diffHeight  = -0.04,
		directional=1,
	)
	
	circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  		xlim = get.cell.meta.data("xlim")
  		ylim = get.cell.meta.data("ylim")
  		sector.name = get.cell.meta.data("sector.index")
  		circos.text(mean(xlim), ylim[1] + .1,
			sector.name,
			facing = "clockwise",
			niceFacing = TRUE,
			adj = c(0, 0.5),
			col = "lightgray")
		},
	bg.border = NA)

	# plot total mutual information content per residue
	totcont = apply(posmat, 2, sum)
	totcont = as.data.frame(cbind(seq(from=1, to=length(totcont)),
					totcont))	
	
	plot(totcont, type='h',
		cex.axis = 2,
		cex.lab = 2,
		cex = 2,
		xlab = "Fragment",
		ylab = "Mutual information content")
	
	contcut = quantile(totcont$totcont, 0.95)
	contrm = apply(totcont, 1, function(row) all(row > contcut ))
	sigcont = totcont[contrm,]


	print(sigcont$V1 + 13)
	
}

pdf("aa_12_circosplt.pdf")
interprotcircosplt(dif, 0.05, 1, 515, 516, 1030, aa12)
dev.off()

pdf("aa_34_circosplt.pdf")
interprotcircosplt(dif, 0.05, 1031, 1545, 1546, 2060, aa34)
dev.off()

pdf("cc_13_circosplt.pdf")
interprotcircosplt(dif, 0.05, 1, 515, 1031, 1545, cc13)
dev.off()

pdf("cc_24_circosplt.pdf")
interprotcircosplt(dif, 0.05, 516, 1030, 1546, 2060, cc24)
dev.off()







