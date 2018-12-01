#==============================================================================
# AlloHubMat
#' Compute the covariance overlap between each of the trajectory blocks from an input list
#' of eigensystems.
#'
#' @param eigensys.list A list of eigen-systems determined from comp_eigensystem()
#' @param eigfrom First eigenvector considered in computing the covariance overlap
#' @param eigto Last eigenvector considered in computing the covariance overlap
#' @return A matrix containing the covariance overlap between each of the trajectory blocks
#' @export 
#==============================================================================

block_overlap = function(eigensys.list, eigfrom = 1, eigto = 10){

	# parse default values if user does not specify
	eigfrom = as.numeric(ifelse(is.na(eigfrom), 1, eigfrom))
	eigto = as.numeric(ifelse(is.na(eigto), 10, eigto))

	# number of blocks in trajectory	
	nBlock = length(eigensys.list)

	# matrix to store overlaps between all block pairs
	traj.overlap = matrix(0, nrow = nBlock, ncol = nBlock)
	
	# range of eigenvalues used in calculating the covariance overlap
	eigrange = eigfrom:eigto


	print("COMPUTING COVARIANCE OVERLAP OF MUTUAL INFORMATION MATRICES")
	
	# measure time taken to execute 
	system.time(

		# compute the covariance overlap between each combination of 
		# eigen-systems

		for (i in 1:(nBlock-1)) {
			for (j in (i+1):nBlock) {
				traj.overlap[i, j] = covar_overlap(eigensys.list[[i]], eigrange, eigensys.list[[j]], eigrange)

				# the overlap between matrix i and matrix j is identical to the overlap between matrix j and matrix i
				traj.overlap[j, i] = traj.overlap[i, j];
			}

		}

	)

	# return matrix containing covarince overlap values
	return(traj.overlap)
}

#==============================================================================
