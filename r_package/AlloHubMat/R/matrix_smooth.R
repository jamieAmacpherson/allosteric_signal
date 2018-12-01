#===============================================================================
# AlloHubMat
#' Smooth the matrix of covariance overlaps using a kernal smoothing algorithm.
#' 
#' This function takes an input matrix and performs a kernal smoothing.
#' The input is a single
#' matrix and returns a single matrix with the same dimensions as that of the input matrix.
#'
#' @param overlap.matrix Matrix of covariance overlaps, computed using block_overlap()
#' @return A matrix containing the covariance overlap between each of the trajectory blocks
#' @export
#===============================================================================

matrix_smooth = function(overlap.matrix){

	# set diagonal to '0' to limit overlap range to observed values
	diag(overlap.matrix) = 0

	# vectorise matrix values
	overlap.matrix.v = as.vector(overlap.matrix)

	## matrix indices of vector values
	nx = dim(overlap.matrix)[1]
	ny = dim(overlap.matrix)[2]
	x = rep(1:nx, ny)
	y = rep(1:ny, each = nx)
	grid.xy = as.data.frame(cbind(x,y))
	
	## smooth vector values given grid indices, resulting in a list
	overlap.matrix.s = fields::smooth.2d(overlap.matrix.v,
		ind = grid.xy,
		nrow = nx,
		ncol = ny,
		theta = 2)

	# return smoothed matrix
	return(overlap.matrix.s)
}

#===============================================================================

