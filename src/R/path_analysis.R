#! /usr/bin/R

#===============================================================================
# Extract the most significant fragments from a mutual information matrix
# and analyse them using graph theory
# This CUDA version uses the gputools for accelerated matrix operations.
#===============================================================================
## Conformational Sampling and Dynamics of Membrane Proteins
##   From 10-Nanosecond Computer Simulations,
##   Faraldo-Gómez, Sansom et al.; DOI: 10.1002/prot.20257.
## see equations (2) and (3) therein
## It turns out that the normalised function 'omega' yields a more diverse range
##   of overlap values and therefore seems to be more suitable for splitting
##   trajectories into ergodic sectors (which appear as blocks with high
##   subspace overlap).
#______________________________________________________________________________

args = commandArgs(TRUE);
print("Usage: Rscript volcano_plt.R <average matrix> <output prefix>");
print("<average matrix> : Mutual information matrix");
print("<output prefix> : output prefix");

#______________________________________________________________________________
# Load packages
library(reshape2)
library(ggplot2)
library(ggrepel)
library(igraph)
library(ggpubr)
library(bio3d)
library(raster)
library(scales)
library(stringr)
library(foreach)
library(doMC)

## Register the number of CPUs
registerDoMC(6)

#______________________________________________________________________________
# split the matrix into a list of four equal size matrices for each chain


splitmat = function(nmi.matrix.dat, outprefix){
	# read mutual information matrix
	nmi.matrix = as.matrix(read.table(nmi.matrix.dat))
	
	# rename rows and columns
        row.names(nmi.matrix) = c(sprintf("A%s",seq(1:515)),
                                sprintf("B%s",seq(516:1030)),
                                sprintf("C%s",seq(1031:1545)),
                                sprintf("D%s",seq(1546:2060)))

        colnames(nmi.matrix) = c(sprintf("A%s",seq(1:515)),
                                sprintf("B%s",seq(516:1030)),
                                sprintf("C%s",seq(1031:1545)),
                                sprintf("D%s",seq(1546:2060)))

        # extract the upper triangle of the MI matrix
    #nmi.matrix[lower.tri(nmi.matrix)] = 0
	#diag(nmi.matrix) = 0
	
	# generate a list of matrices for each chain
	chainA = nmi.matrix[1:515, 1:515]
	chainB = nmi.matrix[516:1030, 516:1030]
	chainC = nmi.matrix[1031:1545, 1031:1545]
	chainD = nmi.matrix[1546:2060, 1546:2060]

	savematrix = function(matrixdat, chainID){
		write.table(matrixdat,
				file = paste(paste(outprefix, chainID, '.dat', sep='')),
				col.names=F,
				row.names=F)
}
	savematrix(chainA, 'chainA')
	savematrix(chainB, 'chainB')
	savematrix(chainC, 'chainC')
	savematrix(chainD, 'chainD')


	listchains = list(chainA, chainB, chainC, chainD)

	return(listchains)
}

#______________________________________________________________________________
## Process the list of matrices and plot the chain-averaged mutual information
# as a volcano plot 
mat.process = function(listmatrices, sign.pval, sig.nMIval, outprefix, ...){
	#______________________________________________________________________________
	## Arguments:
	# 1. List of mutual information matrices, split into chains
	# 2. P value indicating a significance threshold (eg. p = 0.05)
	# 3. Mutual information significance threshold
	# 4. Output prefix
	#______________________________________________________________________________

	# compute the element-wise mean
	listmean = apply(simplify2array(listmatrices), 1:2, mean)
	listsd = apply(simplify2array(listmatrices), 1:2, sd)


	## Compute the element-wise p-value using a Kruskall Wallis test for significance
	## Function is written to accomodate a list of four matrices
	kwtest = function(input){
		# Input is element-wise matrices
		g = factor(rep(1:4,
			c(1,1,1,1)),
			labels=c('mat1', 'mat2', 'mat3', 'mat4'))

	kw.out = kruskal.test(input, g)$p.value

	return(kw.out)
	}

	# parse the element-wise entries of the list of matrices into the kwtest() function
	# to compute significance 
	#listpval = -log10(apply(simplify2array(listmatrices), 1:2, function(x) {t.test(x)$p.value}))



	is.nan.data.frame <- function(x) {do.call(cbind, lapply(x, is.nan))}
	listsd[is.nan(listsd)] <- 0

	# define the fragment names
	fragnames = paste(melt(listmean)$Var1, melt(listmean)$Var2, sep='_')

	# combine p value and mean matrices
	pltdat = as.data.frame(cbind(melt(listmean)$value,
				(melt(listsd)$value)))

	names(pltdat) = c('MI', 'pval')	

	## plot the result
	# set a significance threshold
	sig.thresh = -log10(sign.pval)
	pav.thresh = sig.nMIval
	nav.thresh = -sig.nMIval 

	# identify points which have a significant change in the mutual information and that
	# that change in mutual information is signficant over the four chains	
	pltdat$threshold = as.factor((pltdat$MI > pav.thresh | pltdat$MI < nav.thresh) & pltdat$pval > sig.thresh) 

	pltdat$fragnames = paste(melt(listmean)$Var1, melt(listmean)$Var2, sep='_')

	plt = ggplot(pltdat, aes(MI, pval)) +
		geom_point(alpha = 0.5, aes(colour = threshold), size = 0.2) +
		scale_color_manual(values = c('grey', 'forestgreen')) +
		theme_bw() +
	#	xlim(0, 0.4) +
		xlab(expression(paste(Delta^{'FBP'}, ' nMI'))) +
		ylab(expression(paste(-log[10], ' p-value'))) +
		theme(plot.margin = (unit(c(0.2, 0.2, 0.5, 0.5), 'lines')),
			axis.text = element_text(size = rel(2), colour = 'black'),
			axis.title = element_text(size = rel(2), colour = 'black'),
			legend.position = 'none') +

		geom_text_repel(data = subset(pltdat, (pltdat$MI > 0.2 | pltdat$MI < -0.2 ) & pltdat$pval > sig.thresh),
			aes(label = fragnames),
			size = 2,
			colour = 'black',
			box.padding = unit(0.5, 'lines'),
			point.padding = unit(0.5, 'lines'))

	pdf(paste(outprefix, '_volcano_plt.pdf', sep=''), height = 4, width = 4)
	print(plt)
	dev.off()
}

#mat.process(splitmat(args[1]), 0.07, 0.1, args[2])

#______________________________________________________________________________
## Statistcal plots of mutual information
# plot mean nMI vs. variance of nMI
plt_sig_mu_MI = function(listmatrices, outprefix, ...){

		# compute the element-wise mean
	listmean = apply(simplify2array(listmatrices), 1:2, mean)
	listsd = apply(simplify2array(listmatrices), 1:2, sd)

	is.nan.data.frame <- function(x) {do.call(cbind, lapply(x, is.nan))}
	listsd[is.nan(listsd)] <- 0

	# define the fragment names
	fragnames = paste(melt(listmean)$Var1, melt(listmean)$Var2, sep='_')

	# combine p value and mean matrices
	pltdat = as.data.frame(cbind(melt(listmean)$value,
				(melt(listsd)$value)))

	names(pltdat) = c('MI', 'pval')	

	plt = ggplot(pltdat, aes(MI, pval)) +
		geom_point(alpha = 0.5, size = 0.2) +
		theme_bw() +
		xlab(expression(paste(mu, Delta^{'FBP'}, ' nMI'))) +
		ylab(expression(paste(sigma, Delta^{'FBP'}, ' nMI'))) +
		theme(plot.margin = (unit(c(0.2, 0.2, 0.5, 0.5), 'lines')),
			axis.text = element_text(size = rel(2), colour = 'black'),
			axis.title = element_text(size = rel(2), colour = 'black'),
			legend.position = 'none')

	pdf(paste(outprefix, '_sig_mu.pdf', sep=''), height = 4, width = 4)
	print(plt)
	dev.off()

}
#plt_sig_mu_MI(splitmat('diff_mat.dat'), 'diff')
#plt_sig_mu_MI(splitmat('fbp_nMI.mat', 'fbp'), 'fbp')
#plt_sig_mu_MI(splitmat('apo_nMI.mat', 'apo'), 'apo')


# plot mean nMI vs. coefficient of variance of nMI
plt_cv_mu_MI = function(listmatrices, outprefix, ...){

	# compute the element-wise mean
	listmean = apply(simplify2array(listmatrices), 1:2, mean)
	listsd = apply(simplify2array(listmatrices), 1:2, cv)**2

	is.nan.data.frame <- function(x) {do.call(cbind, lapply(x, is.nan))}
	listsd[is.nan(listsd)] <- 0

	# define the fragment names
	fragnames = paste(melt(listmean)$Var1, melt(listmean)$Var2, sep='_')

	# combine p value and mean matrices
	pltdat = as.data.frame(cbind(melt(listmean)$value,
				(melt(listsd)$value)))

	names(pltdat) = c('MI', 'pval')	

	plt = ggplot(pltdat, aes(MI, pval)) +
		geom_point(alpha = 0.5, size = 0.2) +
		theme_bw() +
		xlab(expression(paste(mu, Delta^{'FBP'}, ' nMI'))) +
		ylab(expression(paste('cv'^{2}, Delta^{'FBP'}), ' nMI')) +
		theme(plot.margin = (unit(c(0.2, 0.2, 0.5, 0.5), 'lines')),
			axis.text = element_text(size = rel(2), colour = 'black'),
			axis.title = element_text(size = rel(2), colour = 'black'),
			legend.position = 'none') +

		# double log scale
		scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x))) +

		#scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
         #     labels = trans_format("log10", math_format(10^.x))) +

		annotation_logticks()

	pdf(paste(outprefix, '_cv_mu.pdf', sep=''), height = 4, width = 4)
	print(plt)
	dev.off()

}
#plt_cv_mu_MI(splitmat('diff_mat.dat', 'diff'), 'diff')
#plt_cv_mu_MI(splitmat('fbp_nMI.mat', 'fbp'), 'fbp')
#plt_cv_mu_MI(splitmat('apo_nMI.mat', 'apo'), 'apo')

# combine the two conditions (Apo and FBP) and plot nMI vs. coefficient of variance
# of the combined data set

plt_comb_cv = function(apolistmatrix, fbplistmatrix, outprefix){

	
	# melt matrix data into a useable data.frame
	melt_data_vals = function(datin){

		mivals = melt(datin)$value
		return(mivals)
	}


	# melt matrix data names into a useable data.frame
	melt_data_names = function(datin){

		fragnames = paste(melt(datin)$Var1,
				melt(datin)$Var2, sep='_')

		return(fragnames)
	}


	# try T-test, if the data is the same, return NA	
	t.test.p.value = function(...) {

		obj = try(t.test(...), silent=TRUE)

		if (is(obj, "try-error")){
			return(NA)}
		else{
			return(obj$p.value)
		}

	}


	# calculate the coef. of variance from two lists of data
	calcstats = function(fbplistdat, apolistdat){


		fbpdatframe = do.call(cbind, lapply(fbplistdat, melt_data_vals))
		apodatframe = do.call(cbind, lapply(apolistdat, melt_data_vals))

		cvsqrd = c()

		#pval.loop.out = function(col1, col2, col3, col4, col5, col6, col7, col8){
	#		pvaldat = -log10(t.test.p.value(
#				c(col1, col2, col3, col4),
#				c(col5, col6, col7, col8)))
#
#			return(pvaldat)
#		}

#		cvsqrd = 

		pval.loop.out = foreach(i=1:nrow(fbpdatframe), .combine=rbind) %dopar% {
		#for (i in c(1:nrow(fbpdatframe))){
			print(paste(i, nrow(fbpdatframe), sep='/'))

			cvdat = -log10(t.test.p.value(apodatframe[i,], fbpdatframe[i,]))

			cvsqrd = append(cvsqrd, cvdat)
		}

		cvsqrd = pval.loop.out
		print(head(cvsqrd))

		print('finished pval calculation')

		log2fc = log2(apply(fbpdatframe, 1, mean) / apply(apodatframe, 1, mean))

		# concatenate the data frames together 
		datframe = as.data.frame(cbind(fbpdatframe, apodatframe))

		# compute the mean and the coef. of variance**2
		datmean = apply(datframe, 1, mean)

		# remove NAs from coef. of variance**2 calculation
		cvsqrd[is.na(cvsqrd)] = 0

		# cbind the resulting statistics
		outdat = as.data.frame(cbind(datmean, cvsqrd, log2fc))

		return(outdat)
	}

	dat = calcstats(fbplistmatrix, apolistmatrix)

 	# combine the two conditions into a single data.frame
 	combdatframe = as.data.frame(cbind(melt_data_names(
 		fbplistmatrix[1]),
 		dat))

 	print(head(combdatframe))

 	names(combdatframe) = c('fragnames', 'mu', 'pval', 'log2fc')

	# set a significance threshold
	sig.thresh = -log10(0.01)
	pav.thresh = 2
	nav.thresh = -2

	# identify points which have a significant change in the mutual information and that
	# that change in mutual information is signficant over the four chains	
	combdatframe$threshold = as.factor((combdatframe$log2fc > pav.thresh | combdatframe$log2fc < nav.thresh) & combdatframe$pval > sig.thresh) 
 	

	## plot mutual information mean vs. coefficient of variance
 	
 	plt = ggplot(combdatframe, aes(x = log2fc, y = pval)) +
		
		# points for FBP
		geom_point(size = 0.1, aes(colour = threshold)) +
		scale_color_manual(values = c('grey', 'forestgreen')) +
		
		# theme settings
		theme_bw() +
		xlab(expression(paste(log[2], ' fold-change'))) +
		ylab(expression(paste(-log[10], ' p-value'))) +
		theme(plot.margin = (unit(c(0.2, 0.2, 0.5, 0.5), 'lines')),
			axis.text = element_text(size = rel(2), colour = 'black'),
			axis.title = element_text(size = rel(2), colour = 'black'),
			legend.position = 'none')


		#geom_text_repel(data = subset(combdatframe, (combdatframe$log2fc > pav.thresh) & combdatframe$pval > sig.thresh),
		#	aes(label = fragnames),
		#	size = 2,
		#	colour = 'black',
		#	box.padding = unit(0.5, 'lines'),
		#	point.padding = unit(0.5, 'lines'))


	#pdf(paste(outprefix, '_volcano.pdf', sep=''), height = 4, width = 4)
	#print(plt)
	#dev.off()

	#______________________________________________________________________________
	## positive signficant fragments
	#______________________________________________________________________________
	pos_sign = subset(combdatframe, (combdatframe$log2fc > pav.thresh) & combdatframe$pval > sig.thresh)
	neg_sign = subset(combdatframe, (combdatframe$log2fc < nav.thresh) & combdatframe$pval > sig.thresh)
	

	# order the signficant fragments according to their mutual information value
	extract.sig.frags = function(subsetdat, outstring){

		sig_frags = subsetdat[order(subsetdat$mu),]$fragnames

		ordered_frags.i = c()
		ordered_frags.j = c()
		
		for (i in c(1:(length(sig_frags)-1))){
			ordered_frags.i = append(ordered_frags.i, as.numeric(gsub('[^0-9]', '', str_sub(sig_frags[i],1,4))))

		# depending on the length of the character fragment coupling, extract
		# the correct number of values to give the second fragment in the coupled pair
		if (nchar(as.character(sig_frags[i])) == 9){
			ordered_frags.j = append(ordered_frags.j, as.numeric(gsub('[^0-9]', '', str_sub(sig_frags[i],-5))))
		}

		if (nchar(as.character(sig_frags[i])) == 8){
			ordered_frags.j = append(ordered_frags.j, as.numeric(gsub('[^0-9]', '', str_sub(sig_frags[i],-4))))
		}

		if (nchar(as.character(sig_frags[i])) == 7){
			ordered_frags.j = append(ordered_frags.j, as.numeric(gsub('[^0-9]', '', str_sub(sig_frags[i],-3))))
		}
	}

	# combine the coupled fragment pair and remove duplicate fragments
	ordered_frags = as.data.frame(cbind(ordered_frags.i, ordered_frags.j)) #+ 12
	
	
	# add 12 because structure starts at 12th residue
	#ordered_frags = unique(as.vector(t(ordered_frags))) + 12
	

	# extract the top ten fragments
	top10.frags = tail(ordered_frags, 20)

	# write pymol script to highlight the top 10 identified fragments
	write.table(top10.frags,
		file = paste(outstring, '_hubs_script.txt', sep=''),
		col.names = FALSE,
		row.names = FALSE,
		sep = '+')

	return(list(ordered_frags, top10.frags))

	}
	
	pos_hubs = extract.sig.frags(pos_sign, 'positive')

	neg_hubs = extract.sig.frags(neg_sign, 'negative')

	return(pos_hubs)
	
}

combdat = plt_comb_cv(splitmat('fbp_nMI.mat', 'fbp'), splitmat('apo_nMI.mat', 'apo'), 'catdat')

plt_comb_cv(listx, listy, 'test')


#______________________________________________________________________________
# Calculate the evolutionary conservation for each of the amino acid residues
# contained within each of the hub fragments
#
#



#______________________________________________________________________________
## Calculate the distance for the shortest path between two vertices using the
# Dijkstra algorithm. This is applied over all combinations of vertices between
# the FBP binding pocket and the active site. Residues at the two pockets are 
# defined from Dombrauckas J, et al. Biochemistry 2005, 44, 9417-9429. 
dijsktra.distance = function(matrixdat, pdb.file, ...){
	#______________________________________________________________________________
	## Arguments:
	# 1. Single mutual information matrix (either Apo or +FBP) 
	# 2. pdb file   
	# 3. Extras ...
	#______________________________________________________________________________
	#
	## scale the mutual information matrix by the C-alpha distance matrix
	#
	# compute the C-alpha distance matrix
	distance.mat = function(pdbdat){
		# read pdb file
		print('READING PDBFILE')
		pdb = read.pdb(pdbdat)
		# determine the c-alpha distance matrix
		print('COMPUTING DISTANCE MATRIX')
		distmat = dm(pdb, inds='calpha', mask.lower=FALSE)
		# subset the distance matrix for the chain
		distmat.chain = distmat[1:515, 1:515]
		# return the distance matrix for the chain
		return(distmat.chain)
	}

	# Compute the distance matrix
	distmat = distance.mat(pdb.file)

	## To pre-process the mutual information matrix, the 
	## matrix is first inflected about the x-axis as so:
	## f(x) = 1 - x. This is so that the 'strong' couplings
	## have small values approaching 0 and the 'weak'
	## couplings have values approaching 1. This inflected
	## matrix is then scaled by the distance matrix of the 
	## protein.

	# transpose and translate the mutual information matrix
	# so that values approaching 0 are strong correlations
	# and ~1 are insignificant correlations
	matrixdat.T = 1 - matrixdat
	
	# scale the mutual information matrix by the c-alpha 
	# distance matrix
	scaled.MImatrix = matrixdat.T * distmat

	# remove the lower triangle and diagonal elements of the 
	# mutual information matrix
	scaled.MImatrix[lower.tri(scaled.MImatrix)] = 0
	diag(scaled.MImatrix) = 0
	
	# vector containing FBP-binding pocket fragments
	# as defined in the following publication:
	# Dombrauckas J, et al. Biochemistry 2005, 44, 9417-9429
	fbp.pocket = c(seq(from=419, to=427),
			seq(from=500, to=510),
			seq(from=468, to=478))

	# vector containing active-site fragments
	catalytic.pocket = c(63, 66, 101, 102, 108, 258, 260, 284, 316, 350)


	# convert matrix into weighted graph
	ig = graph.adjacency(scaled.MImatrix,
			mode = 'undirected',
			weighted = TRUE)

	## loop through possible combinations of FBP binding pocket fragments
	## and active site residues
	
	# generate a table of all unique combinations
	combins = as.matrix(expand.grid(fbp.pocket, catalytic.pocket))

	# empty vector to store the distance results
	path.distances.list = list()

	# loop through all possible combinations and calculate the shortest
	# distance between fragments
	for (i in 1:nrow(combins)){
		
		print(paste(i, nrow(combins), sep='/'))
	
		# calculate the distance and store in temp variable
		dist.tmp = distances(ig,
				v = combins[[i,1]],
				to = combins[[i,2]],
				algorithm = 'dijkstra')
		
		path.distances.list[[i]] = dist.tmp # add it to your list
	}
	
	# combine distances into a dataframe
	path.distances = do.call(rbind, path.distances.list) 

	# return the distances
	return(path.distances)
}



#______________________________________________________________________________
## calculate the shortest distances for each MI matrix
calc.distances = function(matdatname, outprefix, pdb.file){
	# read matrix and split it into chains
	matlist = splitmat(matdatname, outprefix)
	
	# calculate the shortest distances over each of the 
	# matrices in the list
	distance.values = lapply(matlist, dijsktra.distance, pdb.file)
	
	return(distance.values)
}

plt.path.distances = function(fbpmat.file, apomat.file, pdb.file){
	## Plot the path distances for the Apo and FBP matrices
	fbpmatlist = calc.distances(fbpmat.file, 'fbp', pdb.file)
	apomatlist = calc.distances(apomat.file, 'apo', pdb.file)
	subt = function(x){return(x-5)}
	fbpmatlist = lapply(fbpmatlist, subt)

	dfdist = length(apomatlist[[1]][,1])

	## unlist lists of distances
	dat = as.data.frame(rbind(
			cbind(as.vector(apomatlist[[1]][,1]),
					rep('chain1'),
					rep('Apo'),
					rep(1, dfdist)),
			cbind(as.vector(fbpmatlist[[1]][,1]),
					rep('chain1'),
					rep('FBP'),
					rep(2, dfdist)),
			cbind(as.vector(apomatlist[[2]][,1]),
					rep('chain2'),
					rep('Apo'),
					rep(1, dfdist)),
			cbind(as.vector(fbpmatlist[[2]][,1]),
					rep('chain2'),
					rep('FBP'),
					rep(2, dfdist)),
			cbind(as.vector(apomatlist[[3]][,1]),
					rep('chain3'),
					rep('Apo'),
					rep(1, dfdist)),
			cbind(as.vector(fbpmatlist[[3]][,1]),
					rep('chain3'),
					rep('FBP'),
					rep(2, dfdist)),
			cbind(as.vector(apomatlist[[4]][,1]),
					rep('chain4'),
					rep('Apo'),
					rep(1, dfdist)),
			cbind(as.vector(fbpmatlist[[4]][,1]),
					rep('chain4'),
					rep('FBP'),
					rep(2, dfdist))))

	names(dat) = c('distance', 'chain', 'system', 'group')

	dat$distance <- as.numeric(as.character(dat$distance))
	dat$group <- as.numeric(as.character(dat$group))

	# statistical comparisons 
	my_comparisons = list(c('1', '2'))
	
	plt = 	ggboxplot(dat, x='group', y='distance', fill='system', size=0.2,
			palette = c("red", "forestgreen"), width=0.6) +
              	  geom_jitter(position=position_jitter(width=.1, height=0), size=0.1) +
              	  theme_bw() +
                  
                  theme(text = element_text(size=20),
                          axis.text.x = element_text(hjust=1)) +
                  facet_wrap( ~ chain, scales='free') +
                  stat_compare_means(comparisons = my_comparisons, map_signif_level = F) +
                  expand_limits(y=0) +
                  ylim(0, 61)

	pdf('graph_distances.pdf')
	print(plt)
	dev.off()

}


plt.path.distances('fbp_nMI.mat', 'apo_nMI.mat', 'prot_ca.pdb')



#______________________________________________________________________________
## Determine the vertices along the shortest path using the. This is applied
# over all combinations of vertices between the FBP binding pocket and the
# active site. Residues at the two pockets are defined from Dombrauckas J,
# et al. Biochemistry 2005, 44, 9417-9429. 
path.trace = function(matrixdat, pdb.file){
	#______________________________________________________________________________
	# Calculate the shortest path between the allosteric pocket
	# and the active site of PKM2
	#______________________________________________________________________________
	## Arguments:
	# 1. Single mutual information matrix (either Apo or +FBP) 
	# 2. pdb file   
	# 3. Extras ...
	#______________________________________________________________________________
	#
	## scale the mutual information matrix by the C-alpha distance matrix
	#
	# compute the C-alpha distance matrix
	distance.mat = function(pdbdat){
		# read pdb file
		print('READING PDBFILE')
		pdb = read.pdb(pdbdat)
		# determine the c-alpha distance matrix
		print('COMPUTING DISTANCE MATRIX')
		distmat = dm(pdb, inds='calpha', mask.lower=FALSE)
		# subset the distance matrix for the chain
		distmat.chain = distmat[1:515, 1:515]
		# return the distance matrix for the chain
		return(distmat.chain)
	}

	# Compute the distance matrix
	distmat = distance.mat(pdb.file)

	## To pre-process the mutual information matrix, the 
	## matrix is first inflected about the x-axis as so:
	## f(x) = 1 - x. This is so that the 'strong' couplings
	## have small values approaching 0 and the 'weak'
	## couplings have values approaching 1. This inflected
	## matrix is then scaled by the distance matrix of the 
	## protein.

	# transpose and translate the mutual information matrix
	# so that values approaching 0 are strong correlations
	# and ~1 are insignificant correlations
	matrixdat.T = 1 - matrixdat
	
	# scale the mutual information matrix by the c-alpha 
	# distance matrix
	scaled.MImatrix = matrixdat.T * distmat	
	# remove the lower triangle and diagonal elements of the 
	# mutual information matrix
	scaled.MImatrix[lower.tri(scaled.MImatrix)] = 0
	diag(scaled.MImatrix) = 0
	
	# vector containing FBP-binding pocket fragments
	fbp.pocket = c(seq(from=419, to=427),
			seq(from=500, to=510),
			seq(from=468, to=478))

	# vector containing active-site fragments
	catalytic.pocket = c(63, 66, 101, 102, 108, 258, 260, 284, 316, 350)


	# convert matrix into weighted graph
	ig = graph.adjacency(scaled.MImatrix,
			mode = 'undirected',
			weighted = TRUE)

	## loop through possible combinations of FBP binding pocket fragments
	## and active site residues
	
	# generate a table of all unique combinations
	combins = as.matrix(expand.grid(fbp.pocket, catalytic.pocket))

	# empty vector to store the distance results
	path.distances.list = list()

	# loop through all possible combinations and calculate the shortest
	# distance between fragments
	for (i in 1:nrow(combins)){
		
		print(paste(i, nrow(combins), sep='/'))
	
		# calculate the distance and store in temp variable
		dist.tmp = shortest_paths(ig,
				from = combins[[i,1]],
				to = combins[[i,2]],
				weights = NULL,
				output = 'vpath')
		
		path.distances.list[[i]] = dist.tmp # add it to your list
	}
	
	# combine distances into a dataframe
	path.distances = do.call(rbind, path.distances.list) 

	# return the distances
	return(path.distances)
}


#______________________________________________________________________________
## calculate the maximum-sum path for each MI matrix
calc.paths = function(matdatname, outprefix, pdb.file){
	# read matrix and split it into chains
	matlist = splitmat(matdatname, outprefix)
	
	# calculate the shortest distances over each of the 
	# matrices in the list
	paths = lapply(matlist, path.trace, pdb.file)
	
	return(paths)
}

## calculate the maximum sum paths for all FBP chains
fbppaths = calc.paths('fbp_nMI.mat', 'fbp', 'prot_ca.pdb')
apopaths = calc.paths('apo_nMI.mat', 'apo', 'prot_ca.pdb')

#______________________________________________________________________________
## Plot the paths onto a structure using the bio3d plugin
plot.paths = function(pathlist, chain){

	pathlist.chain = pathlist[[chain]]

	# read pdb file
	# for each of the following active-site residues:
	# 63, 66, 101, 102, 108, 258, 260, 284, 316, 350
	as.residues = c(63, 66, 101, 102, 108, 258, 260, 284, 316, 350)

	# initialise a list to take the paths from each of the active 
	# site residues
	as.paths = list()

	## there are 63 path combinations for each of the active site residues
	## loop through the list of paths in increments of 31 and extract
	## the vertices

	# initialise variable

	i = 1

	while (i < (nrow(pathlist.chain)-1) ){
		# start of the list
		start.pos = i
		end.pos = i + 31

		print(start.pos)
		print(end.pos)

		# extract the list of vertices
		dat = as.data.frame(
			table(
				unlist(
					pathlist.chain[start.pos:end.pos]) + 12
				)
			) $ Var1

		# convert the string of vertices from class factor to 
		# class numeric
		dat.num = as.numeric(as.character(dat))

		# add 32 positions to the iterator
		list.pos = end.pos/32

		# add the numeric vector of vertices to a list
		as.paths[[list.pos]] = dat.num

		i = i + 32
	}

	## generate a list of strings to draw bonds in PyMOL
	
	# bond selection command, given two input atoms and a name for the 
	# bond - returned as a string
	pymol.bond.sele.command = function(atom1, atom2, nbond, chainID){

		# chain switches
		if(chainID == 1){
			pymolchain = 'A'
		}
		else if(chainID == 2){
			pymolchain = 'B'
		}
		else if(chainID == 3){
			pymolchain = 'C'
		}
		else if(chainID == 4){
			pymolchain = 'D'
		}

		# return PyMOL command to draw a bond between atoms 1 and 2
		return(paste(paste('distance ', paste('mydist', nbond, ', ', sep=''),
			'chain ', pymolchain, ' and ', atom1, '/CA, ', 'chain ', pymolchain, ' and ', atom2, '/CA', sep='')))
	}

	## remove the labels from the selected distances
	pymol.bond.lab.rm.command = function(nbond, chainID){

		# return PyMOL command to remove the distance label drawn for distance element x
		return(paste('hide labels, mydist', nbond, sep=''))
	}


	## function which writes a pymol selection command for pathway
	pymol.bond = function(vertex.vec){
		# initialise a vector to catch the selection string commands
		pymol.commands = c()

		# loop through the vertex vector calling pymol.bond.sele.command() for 
		# each successive combination of path vertices
		for ( w in c(1:length(vertex.vec)-1)){
			cmd.tmp = pymol.bond.sele.command(vertex.vec[w], vertex.vec[w + 1], w, chain)
			pymol.commands = append(pymol.commands, cmd.tmp)
		}

		# remove the first element of the list (redundant output)
		pymol.commands = pymol.commands[-1]
		return(pymol.commands)
	}

	path.commands = lapply(as.paths, pymol.bond)

	## function which writes a pymol selection command for removing the labels from the pathway lines
	pymol.rm.label = function(vertex.vec){
		# initialise a vector to catch the selection string commands
		pymol.commands = c()

		# loop through the vertex vector calling pymol.bond.sele.command() for 
		# each successive combination of path vertices
		for ( w in c(1:length(vertex.vec)-1)){
			cmd.tmp = pymol.bond.lab.rm.command(w)
			pymol.commands = append(pymol.commands, cmd.tmp)
		}

		# remove the first element of the list (redundant output)
		pymol.commands = pymol.commands[-1]
		return(pymol.commands)
	}

	path.rm.lab = lapply(as.paths, pymol.rm.label)
	#return(as.paths)

	## write the PyMOL script file to the working directory
	fileConn = file(paste('as_paths_chain', chain, '_script.txt', sep=''), 'w')
	lapply(path.commands,
		write,
		fileConn,
		append=TRUE)

	close(fileConn)

	## write PyMOL script that removes the distances from the drawn paths in the first script above
	fileConn = file(paste('as_paths_chain', chain, '_script_rmlab.txt', sep=''), 'w')
	lapply(path.rm.lab,
		write, 
		fileConn,
		append=TRUE)

	close(fileConn)	

	#fileConn = file(paste('as_paths_chain', chain, '_script.txt'))
	#writeLines(path.commands, fileConn)
	#close(fileConn)

	print('PYMOL SCRIPT FILE SAVED TO WD, CONTAINING LIST OF RESIDUES IN ALLOSTERIC PATHWAY')
}

plot.paths.allchains = function(nchains){
	# loop through all chains in the structure
	# and extract the allosteric pathways
	for(i in c(1:nchains)){
		# FBP-bound pathways
		plot.paths(fbppaths, i)

		# Apo PKM2 pathways
		plot.paths(apopaths, i)

	}
}

plot.paths.allchains(4)



identify.hub.resids = function(){
	# identify hub residues along each of the 10 paths
	hub.resids = c(491, 247, 357, 197, 492, 413, 237, 70, 307, 178)

	# active site residues
	as.resids = c(75, 78, 113, 114, 120, 270, 272, 296, 328, 362)

	## for all combinations of hub and path residue
	# generate a table of all unique combinations
	combins = as.matrix(expand.grid(hub.resids, as.resids))

	# empty vector to store the distance results
	path.distances.list = list()

	# loop through all possible combinations and calculate the shortest
	# distance between fragments
	for (i in 1:nrow(combins)){
		
		print(paste(i, nrow(combins), sep='/'))
	
		# calculate the distance and store in temp variable
		dist.tmp = distances(ig,
				v = combins[[i,1]],
				to = combins[[i,2]],
				algorithm = 'dijkstra')
		
		path.distances.list[[i]] = dist.tmp # add it to your list
	}
	
	# combine distances into a dataframe
	path.distances = do.call(rbind, path.distances.list) 

}
	

}
test = plot.paths(fbppaths, 1)

#______________________________________________________________________________
## print maximum sum path of mutual information matrix 

maxsum.path = function(matrix.dat, start.node, end.node, pathlength){	
	#______________________________________________________________________________
	## Arguments:
	# 1. Mutual information matrix, diagonal and lower triangle set equal to zero 
	# 2. Start vertex  
	# 3. End vertex
	#______________________________________________________________________________

	# set lower triangle and diagonal of matrix to zero	
	matrix.dat[lower.tri(matrix.dat)] = 0
	diag(matrix.dat) = 0
	
	#
	# stop function if the diagonal is not zero
	for ( i in diag(matrix.dat)){
		if(i != 0){
			stop(call. = FALSE)
			print('set diagonal of mutual
			information matrix to zero')
		}
	}

	# stop funciton if lower triangle elements of matrix
	# are greater than zero
	for ( i in matrix.dat[lower.tri(matrix.dat)]){
		if(i != 0){
			stop(call.=FALSE)
		}
	}


	# reshape the matrix so that the first elemnt is the start node and the 
	# last element is the end node
	chart = matrix.dat[start.node:end.node, start.node:end.node]	

	# keep the original chart
	original = chart

	# find the maximum sum path between the two nodes
	# preset the pathlength
	n = pathlength


	chart[2, 1:2] <- chart[2, 1:2] + chart[1, 1]
	for (i in 3:n) {
		chart[i, 1] <- chart[i, 1] + chart[(i - 1), 1]
		chart[i, i] <- chart[i, i] + chart[(i - 1), (i - 1)]
  
		for (j in 2:(i - 1)) {
			chart[i, j] <- chart[i, j] + max(chart[(i - 1), (j - 1):j])
		}
	}


	result <- max(chart[n, ])

	#
	# get the path
	#
	# route will have n elements
	route <- rep(0,n) 

	# index of last max (in last row)
	route[n] = which(chart[n,]==max(chart[n,],na.rm=TRUE))[[1]] 
	
	route[1] = 1 # top of the pyramid
	for (i in (n-1):2) { # starting from bottom, going to top
	
		left <- if (route[i+1] > 1) route[i+1] -1 else 1 # index of left element
		right <- if (route[i+1] < i+1) route[i+1] else i # index of right element
		route[i] = which(chart[i,]==max(chart[i,left], chart[i,right], na.rm=TRUE))[[1]] # choose the higher value
	}


	# return the results
	print(paste('result: ', result))
	print(paste('Route: ', route))

	# check the route
	checksum = 0;
	for (i in 1:n) checksum = checksum + original[i, route[i]]
	
	print(paste('Checksum: ', checksum))

}


trace.paths = function(matdatname, outprefix){
        # read matrix and split it into chains
        matlist = splitmat(matdatname, outprefix)

        # vector containing FBP-binding pocket fragments
        fbp.pocket = c(seq(from=419, to=427),
                        seq(from=500, to=510),
                        seq(from=468, to=478))

        # vector containing active-site fragments
        catalytic.pocket = c(63, 66, 101, 102, 108, 258, 260, 284, 316, 350)

        # calculate the shortest distances over each of the 
        # matrices in the list
        distance.values = lapply(matlist, maxsum.path, 419, 66, 200)

        return(distance.values)
}


#______________________________________________________________________________
## Compute the decay of the mutual information signal proximal to allosteric
## hubs

MIdecay = function(MImatrix){
	print('''
		Compute the decay of the mutual information signal proximal to allosteric hubs.

		Inputs:
		(1) Mutual information matrix

		Output:
		(2) Autocorrelation function decay of mutual information signal either side of 
		allosteric fragment hubs.
		''')

	# hard-coded allosteric hubs
	allosteric.hubs = c(109, 420, 191, 293, 292, 312, 343, 477, 232, 274)



}



