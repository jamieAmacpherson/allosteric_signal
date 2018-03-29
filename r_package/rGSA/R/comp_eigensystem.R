#' Compute the eigensystem for the list of mutual information matrices
#'
#' This function computes the covariance matrix and eigensystems for each of
#' the mutual information matrices, provided to the function, in the input list.
#'
#' @param MImatlist A list of mutual information matrices
#' @return A list containing the eigensystems for each of the mutual information matrices
#' in the input list
#' @export 



comp_eigensystem = function(MImatlist){

	# number of blocks in trajectory
	nBlock = length(MImatlist)

	# initialise list for the covariance matrices and eigensystems
	covtraj = list(nBlock)
	eigtraj = list(nBlock)

	# compute the covarince matrix for each of the elements in the input list
	print("COMPUTING COVARIANCE MATRICES")
	covtraj = lapply(MImatlist, covar)

	# conpute eigensystem for each matrix in input list
	print("COMPUTING EIGEN-SYSTEMS")
	eigtraj = lapply(covtraj, eigen)

	return(eigtraj)	
}