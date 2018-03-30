#' Smooth the matrix of covariance overlaps using a kernal smoothing algorithm.
#' 
#' This function takes an input matrix and performs a kernal smoothing. The input is a single
#' matrix and returns a single matrix with the same dimensions as that of the input matrix.
#'
#' @param overlap.matrix Matrix of covariance overlaps, computed using block_overlap()
#' @param eigfrom First eigenvector considered in computing the covariance overlap
#' @param eigto Last eigenvector considered in computing the covariance overlap
#' @return A matrix containing the covariance overlap between each of the trajectory blocks
#' @export
