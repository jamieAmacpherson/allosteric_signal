#===============================================================================
# AlloHubMat
#' Compute the covariance overlap between two matrices
#'
#' This function computes the covariance overlap between two matrices of identical
#' dimensions.
#'
#' The matrix calculation is performed according to the following equation (LaTex):
#' 		\begin{equation}
#' 			\Omega_{A:B} = 1 - \lbrace \frac{\sum_{i=1}^{3N}(\lambda_{i}^{A} + \lambda_{i}^{B}) - 2 \sum_{i=1}^{3N} \sum_{j=1}^{3N} (\lambda_{i}^{A} \lambda_{i}^{B})^{0.5} (\textbf{v}_{i}^{A} \cdot \textbf{v}_{i}^{B})^{2}}{\sum_{i=1}^{3N}(\lambda_{i}^{A} + \lambda_{i}^{B})} \rbrace ^{0.5}
#' 		\end{equation}
#'
#'
#' @param esA First eigen-system
#' @param irange Eigenvalue range of first eigen-system
#' @param esB Second eigen-system
#' @param jrange Eigenvalue range of second eigen-system
#' @return Return the covariance overlap between the two input eigen-systems
#'		   for the defined ranges of eigenvalues
#' @export
#===============================================================================

covar_overlap = function(esA, irange, esB, jrange) {

	## assert index range is within matrix dimension	
	stopifnot(dim(esA$vectors[2]) >= range(irange)[2]);
	stopifnot(dim(esB$vectors[2]) >= range(jrange)[2]);
	
	## assert index ranges are equally long
	stopifnot(length(irange) == length(jrange));

	# sum terms of the equation are called here 'ta' and 'tb' */
	ta = sum(sapply(c(1:length(irange)), function(x) {
		esA$values[irange[x]] + esB$values[jrange[x]]; 
	}));

	tb = sum(sapply(c(1:length(irange)), function(x) {
		sqrt(esA$values[irange[x]] * esB$values[jrange[x]]) *
	  	esA$vectors[ , irange[x]] %*% esB$vectors[ , jrange[x]];		
	}));

	## using the 'abs' function to avoid negatives (from rounding errors)
	outcomp = (1 - sqrt(abs((ta - (2 * tb))) / ta))

	return(outcomp)
}

#===============================================================================

