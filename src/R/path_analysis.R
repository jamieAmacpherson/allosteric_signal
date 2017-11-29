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
        nmi.matrix[lower.tri(nmi.matrix)] = 0
	diag(nmi.matrix) = 0
	
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

	# compute the element-wise p value from non-parametric t test
	listpval = -log10(apply(simplify2array(listmatrices), 1:2, function(x) {t.test(x)$p.value}))

	is.nan.data.frame <- function(x) {do.call(cbind, lapply(x, is.nan))}
	listpval[is.nan(listpval)] <- 0

	# define the fragment names
	fragnames = paste(melt(listmean)$Var1, melt(listmean)$Var2, sep='_')

	# combine p value and mean matrices
	pltdat = as.data.frame(cbind(melt(listmean)$value,
				(melt(listpval)$value)))

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

	while (i < (length(pathlist)-1) ){
		# start of the list
		start = i
		end = i + 31

		# extract the list of vertices
		dat = as.data.frame(
			table(
				unlist(
					pathlist[[chain]][start:end]) + 12
				)
			)

		

	dat = as.data.frame(
		table(
			unlist(
				pathlist[[chain]][begin:end])+12
			)
		)$Var1
	

	print(dat)

	## TODO: plot each path onto the structure:
	# for f in paths:
	#	for vertices in f:
	#		draw line between vertices
	# export as python executable script

	#pymol(pdb,
	#	as = 'putty',
	#	type = 'script',
	#	)

}


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


