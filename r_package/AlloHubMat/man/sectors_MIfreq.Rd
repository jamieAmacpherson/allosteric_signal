% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_sectors.R
\name{sectors_MIfreq}
\alias{sectors_MIfreq}
\title{Plot the frequency distribution of the ergodic sector mutual information.}
\usage{
sectors_MIfreq(sector.list, n.sigmas)
}
\arguments{
\item{sector.list}{A list of matrices, each is an averaged mutual information matrix for an ergodic sector}

\item{n.sigmas}{An integer number of standard deviations used in calculating the mutual information significance threshold (default value of 2)}
}
\value{
A density histogram plot of the mutual information of the ergodic sectors
}
\description{
This function generates a frequency distribution plot of the mutual information contained by the ergodic sectors
identified by extract_sectors(), and calculates a significance cutoff to distinguish significant correlations
from the noise.
}
