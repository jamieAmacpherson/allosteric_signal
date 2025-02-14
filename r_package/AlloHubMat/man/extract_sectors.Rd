% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sectors.R
\name{extract_sectors}
\alias{extract_sectors}
\title{Extract the ergodic sectors from a trajectory of mutual information matrices.}
\usage{
extract_sectors(sector.info, MImatlist)
}
\arguments{
\item{sector.info}{a dataframe containing a summary of any detected ergodic sectors, including the 
time in the MD trajectory in which they occur. Output of detect_sectors().}

\item{MImatlist}{First eigenvector considered in computing the covariance overlap}
}
\value{
A list of matrices, each is an averaged mutual information matrix for an ergodic sector
}
\description{
This function takes the list of mutual information matrices <output of read_mats()> with the detected ergodic
sectors <output of detect_sectors()>, and extracts the ergodic sectors from the trajectory of mutual information
matrices. The function returns a list of averaged ergodic sectors, 
matrix generated by matrix_smooth(), and returns a list of ergodic sectors.
}
