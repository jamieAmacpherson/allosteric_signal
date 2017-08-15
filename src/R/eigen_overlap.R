#! /usr/bin/R

#===============================================================================
# Methods to compute the overlap of eigenspaces and their sub-spaces
#===============================================================================
## Conformational Sampling and Dynamics of Membrane Proteins
##   From 10-Nanosecond Computer Simulations,
##   Faraldo-Gómez, Sansom et al.; DOI: 10.1002/prot.20257.
## see equations (2) and (3) therein

library("bio3d");

#______________________________________________________________________________
## dummy data matrices mat.A and mat.B
mat.A = runif(100);
dim(mat.A) = c(10,10);
## covariance matrix
cov.A = cov(mat.A);
## eigensystem
eig.A = eigen(cov.A);

mat.B = runif(100);
dim(mat.B) = c(10,10);
## covariance matrix
cov.B = cov(mat.B);
## eigensystem
eig.B = eigen(cov.B);

## range of vector indices forming vector spaces in eig.A (and eig.B)
## no overlap
irange.0 = c(1:5); 
jrange.0 = c(6:10);
## complete overlap
irange1 = c(1:5);
jrange.1 = c(1:5);

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
psiAB(eig.A, irange.0, eig.A, jrange.0);
## should be close to 1
psiAB(eig.A, irange1, eig.A, jrange.1);

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
omegaAB(eig.A, irange.0, eig.A, jrange.0);
## should be close to 1
omegaAB(eig.A, irange1, eig.A, jrange.1);

#______________________________________________________________________________
## test data
pdb = read.pdb("./md1_ca.pdb");
print(pdb);
print(pdb$xyz);

mat.dcd = read.dcd("./md1_ca.dcd");
print(mat.dcd);

## covariance matrix
cov.dcd = cov(mat.dcd);
## eigensystem
eig.dcd = eigen(cov.dcd);

## no overlap
irange.0 = c(1:15); 
jrange.0 = c(16:30);
## complete overlap
irange.1 = c(1:30);
jrange.1 = c(1:30);

## subspace overlap 'psi'
psiAB(eig.dcd, irange.0, eig.dcd, jrange.0);
psiAB(eig.dcd, irange.1, eig.dcd, jrange.1);

## covariance overlap 'omega'
omegaAB(eig.dcd, irange.0, eig.dcd, jrange.0);
omegaAB(eig.dcd, irange.1, eig.dcd, jrange.1);


#______________________________________________________________________________
#______________________________________________________________________________
## decapeptide
## PDB structure
pdb = read.pdb("../../test_sys/deca_A_md/vacuum/md/md1_ca.pdb");
print(pdb);
## DCD trajectory
dcd = read.dcd("../../test_sys/deca_A_md/vacuum/md/md1_rottrans_ca.dcd");
print(dcd);

#______________________________________________________________________________
## trajectory is 10ns, 5001 conformers
covtraj = list(100);
eigtraj = list(100);

## create 100 blocks (of 50 conformers) and compute covariance matrices
##   and derive eigensystems
#for (i in 1:100) {
#	inc = (i - 1) * 50;
#	covtraj[[i]] = cov(dcd[(1 + inc):(50 + inc), ]);
#	eigtraj[[i]] = eigen(covtraj[[i]]);
#}

## create  blocks of increasing size (by 50), like 1:50, 1:100, 1:150, ...
for (i in 1:100) {
	covtraj[[i]] = cov(dcd[(1:(i*50)), ]);
	eigtraj[[i]] = eigen(covtraj[[i]]);
}

## compute overlaps for successive block sizes
iranged = c(1:2);
jranged = c(1:2);

psitraj = vector(length = 99);
omegatraj = vector(length = 99);
for (i in 1:99) {
	psitraj[[i]] = psiAB(eigtraj[[i]], iranged, eigtraj[[i+1]], jranged);
	omegatraj[[i]] = omegaAB(eigtraj[[i]], iranged, eigtraj[[i+1]], jranged);
}

plot(psitraj, type = "line");
plot(omegatraj, type = "line");


#===============================================================================

