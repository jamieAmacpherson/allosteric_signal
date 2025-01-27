% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/path_analysis.R
\name{path.trace}
\alias{path.trace}
\title{Determine the allosteric pathway between a user-defined ligand binding fragment and an active-site fragment.}
\usage{
path.trace(sector.list, frag.from, frag.to, pdb.file)
}
\arguments{
\item{sector.list}{A list of matrices, each is an averaged mutual information matrix for an ergodic sector}

\item{frag.from}{Starting fragment in the pathway (a fragment at the allosteric pocket)}

\item{frag.to}{Ending fragment in the pathway (a fragment at the active site)}

\item{pdb}{Topology of the protein in .pdb format.}
}
\value{
A vector containing the fragments along the predicted allosteric pathway
}
\description{
This function determines an allosteric pathway between an user-defined ligand binding and active site fragment.
The function takes the ergodic sector mutual information matrices from extract_sectors(), a fragment at the 
allosteric site and a fragment at the active site, and returns a list containing the fragments connecting a 
minimal-distance pathway.
}
