#! /usr/bin/R

#===============================================================================
# Methods to compute the overlap of eigenspaces and their sub-spaces
#===============================================================================
## Conformational Sampling and Dynamics of Membrane Proteins From 10-Nanosecond Computer Simulations
##  Faraldo-Gómez, Sansom et al.; DOI: 10.1002/prot.20257
##  see equations (2) and (3)

## dummy matrices
matA = runif(100);
dim(matA) = c(10,10);
covA  = cov(matA);
eigA = eigen(covA);

matB = runif(100);
dim(matB) = c(10,10);
covB  = cov(matB);
eigB = eigen(covB);

## range of vector indices forming vector spaces in eigA (and eigB)
## no overlap
irange0 = c(1:5); 
jrange0 = c(6:10);
## complete overlap
irange1 = c(1:5);
jrange1 = c(1:5);

#______________________________________________________________________________
## subspace overlap 'psi' of eigensystem vector spaces A and B formed by vectors
##   with index ranges 'irange' and 'jrange';  equation (2)
psiAB = function(esA, irange, esB, jrange) {
	## assert index range is within matrix dimension	
	stopifnot(dim(esA$vectors[2]) >= range(irange)[2]);
	stopifnot(dim(esB$vectors[2]) >= range(jrange)[2]);
	## assert index ranges are equally long
	stopifnot(length(irange) == length(jrange));

	psiSum = sum(sapply(c(1:length(irange)), function(x) {
	  	esA$vectors[ , irange[x]] %*% esB$vectors[ , jrange[x]];
	  }))
  return(psiSum / length(irange));
}

## should be close to 0
psiAB(eigA, irange0, eigA, jrange0);
## should be close to 1
psiAB(eigA, irange1, eigA, jrange1);

#______________________________________________________________________________
## covariance overlap 'omega'; equation (3)
omegaAB = function(esA, irange, esB, jrange) {
	## assert index range is within matrix dimension	
	stopifnot(dim(esA$vectors[2]) >= range(irange)[2]);
	stopifnot(dim(esB$vectors[2]) >= range(jrange)[2]);
	## assert index ranges are equally long
	stopifnot(length(irange) == length(jrange));

	# sum terms of the equation are called here 'ta' and 'tb' */
	ta = sum(sapply(c(1:length(irange)), function(x) {
		esA$values[irange[x]] + esB$values[jrange[x]]; 
	}))

	tb = sum(sapply(c(1:length(irange)), function(x) {
		sqrt(esA$values[irange[x]] * esB$values[jrange[x]]) *
	  	esA$vectors[ , irange[x]] %*% esB$vectors[ , jrange[x]];		
	}))

	## using the 'abs' function to avoid negatives (from rounding errors)
	return (1 - sqrt(abs((ta - (2 * tb))) / ta));
}

## should be close to 0
omegaAB(eigA, irange0, eigA, jrange0);
## should be close to 1
omegaAB(eigA, irange1, eigA, jrange1);


#===============================================================================

