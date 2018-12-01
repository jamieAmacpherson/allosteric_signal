#===============================================================================
# AlloHubMat
#' Determine the allosteric pathway between a user-defined ligand binding fragment and an active-site fragment.
#' 
#' This function determines an allosteric pathway between an user-defined ligand binding and active site fragment.
#' The function takes the ergodic sector mutual information matrices from extract_sectors(), a fragment at the 
#' allosteric site and a fragment at the active site, and returns a list containing the fragments connecting a 
#' minimal-distance pathway. 
#'
#' @param sector.list A list of matrices, each is an averaged mutual information matrix for an ergodic sector
#' @param frag.from Starting fragment in the pathway (a fragment at the allosteric pocket)
#' @param frag.to Ending fragment in the pathway (a fragment at the active site)
#' @param pdb Topology of the protein in .pdb format.
#' @return A vector containing the fragments along the predicted allosteric pathway
#' @export
#===============================================================================

path.trace = function(sector.list, frag.from, frag.to, pdb.file){


	# if there are more than a single ergodic sector in the list, then average over the ergodic
	# sector mutual information matrices
	if(length(sector.list) > 1){
		
		# perform an element-wise averaging over the list of sector matrices
		# to determine the average mutual information matrix
		mat = Reduce("+", sector.list) / length(sector.list)


	} else{

		# if a single ergodic sector is identified, extract the first (and only) element of the sector.list
		# list as the matrix
		mat = sector.list[[1]]
	}


	# ensure that the user-defined fragments are within the size of the mutual information matrix
	if(dim(mat)[1] < frag.from | dim(mat)[1] < frag.to){

		stop('FATAL ERROR: allosteric fragment and/or active-site fragment is outside of the mutual information matrix')
	}


	## scale the mutual information matrix by the C-alpha distance matrix
	#
	# compute the C-alpha distance matrix
	distance.mat = function(pdbdat){
		# read pdb file
		print('READING PDBFILE')
		pdb = bio3d::read.pdb(pdbdat)

		# determine the c-alpha distance matrix
		print('COMPUTING DISTANCE MATRIX')
		distmat = bio3d::dm(pdb, inds='calpha', mask.lower=FALSE)

		# if the pdb file is not of the dimensions (n+3, n+3) return an error
		# and exit the function
		if(isFALSE( dim(distmat)[1] + 3 == dim(mat)[1] & dim(distmat)[2] + 3 == dim(mat)[2] )){

			stop(
				writeLines('FATAL ERROR: mismatch between the topology (.pdb) file and the dimensions of the mutual information matrix. \n 
					Please ensure that you are using the correct topology file.'))

		}

		# subset the distance matrix for the chain
		distmat.chain = distmat[1:515, 1:515]

		# return the distance matrix for the chain
		return(distmat.chain)
	}

	distmat = distance.mat(pdb.file)

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
	matrixdat.T = 1 - mat
	
	# scale the mutual information matrix by the c-alpha 
	# distance matrix
	scaled.MImatrix = matrixdat.T * distmat	

	# remove the lower triangle and diagonal elements of the 
	# mutual information matrix
	scaled.MImatrix[lower.tri(scaled.MImatrix)] = 0
	diag(scaled.MImatrix) = 0

	# convert matrix into weighted graph
	ig = igraph::graph.adjacency(scaled.MImatrix,
			mode = 'undirected',
			weighted = TRUE)


	# calculate the distance and store in temp variable
	pathway = igraph::shortest.paths(ig,
			from = frag.from,
			to = frag.to,
			weights = NULL,
			output = 'vpath')


	# return the predicted allosteric pathway
	return(pathway)

}

#===============================================================================

