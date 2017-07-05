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

#______________________________________________________________________________
## some toy data: 3 Hydrogen atoms (H1, H2, H3) with 3 trajectories:
##   H1 moves along the room diagonal of Catesian x,y,z
##   H2 does the same, but returns after half of the distance
##   H3 flies along x, then along y, then along z
## trajectory length: 120 steps

H1.x = seq(from = 1, to = 120, by = 1);
H1.y = seq(from = 1, to = 120, by = 1);
H1.z = seq(from = 1, to = 120, by = 1);
H1 = cbind(H1.x, H1.y, H1.z);

H2.x = c(seq(from = 1, to = 60, by = 1), seq(from = 60, to = 1, by = -1));
H2.y = c(seq(from = 1, to = 60, by = 1), seq(from = 60, to = 1, by = -1));
H2.z = c(seq(from = 1, to = 60, by = 1), seq(from = 60, to = 1, by = -1));
H2 = cbind(H2.x, H2.y, H2.z);

H3.x = c(seq(from = 1, to = 40, by = 1), rep(40, 80));
H3.y = c(rep(1, 40), seq(from = 1, to = 40, by = 1), rep(40, 40));
H3.z = c(rep(1, 80), seq(from = 1, to = 40, by = 1));
H3 = cbind(H3.x, H3.y, H3.z);

#______________________________________________________________________________
## total trajectory
H = cbind(H1, H2, H3);

## trajectory split into 2 parts 
H_t2_1 = H[1:40, ];
covH_t2_1 = cov(H_t2_1);
eigH_t2_1 = eigen(covH_t2_1);

H_t2_2 = H[41:80, ];
covH_t2_2 = cov(H_t2_2);
eigH_t2_2 = eigen(covH_t2_2);

## trajectory split into 3 parts
H_t3_1 = H[1:40, ];
covH_t3_1 = cov(H_t3_1);
eigH_t3_1 = eigen(covH_t3_1);

H_t3_2 = H[41:80, ];
covH_t3_2 = cov(H_t3_2);
eigH_t3_2 = eigen(covH_t3_2);

H_t3_3 = H[81:120, ];
covH_t3_3 = cov(H_t3_3);
eigH_t3_3 = eigen(covH_t3_3);

#______________________________________________________________________________
iranget = c(1:2);
jranget = c(1:2);

## t2_1 t2_2
psiAB(eigH_t2_1, iranget, eigH_t2_2, jranget);
omegaAB(eigH_t2_1, iranget, eigH_t2_2, jranget);

## t3_1 t3_2
psiAB(eigH_t3_1, iranget, eigH_t3_2, jranget);
omegaAB(eigH_t3_1, iranget, eigH_t3_2, jranget);

## t3_1 t3_3
psiAB(eigH_t3_1, iranget, eigH_t3_3, jranget);
omegaAB(eigH_t3_1, iranget, eigH_t3_3, jranget);

## t3_2 t3_3
psiAB(eigH_t3_2, iranget, eigH_t3_3, jranget);
omegaAB(eigH_t3_2, iranget, eigH_t3_3, jranget);


#===============================================================================

