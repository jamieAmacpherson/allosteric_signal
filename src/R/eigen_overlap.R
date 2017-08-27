#! /usr/bin/R

#===============================================================================
# Compute the overlap of eigenspaces and their sub-spaces
#===============================================================================
## Conformational Sampling and Dynamics of Membrane Proteins
##   From 10-Nanosecond Computer Simulations,
##   Faraldo-Gómez, Sansom et al.; DOI: 10.1002/prot.20257.
## see equations (2) and (3) therein
## It turns out that the normalised function 'omega' yields a more diverse range
##   of overlap valuse and therefore seems to be more suitable for dissecting
##   trajectories into ergodic sectors (which appear as blocks with high
##   subspace overlap).

#______________________________________________________________________________
## LIBRARIES and FUNCTIONS
#______________________________________________________________________________

library("bio3d");

args = commandArgs(TRUE);
print("Usage: Rscript eigen_overlap.R <pdb> <dcd> <blocksize> <eigfrom> <eigto>"); 
print("<pdb> : molecular structure file in PDB format");
print("<dcd> : molecular trajectory file in DCD format");
print("<blocksize> : number of conformers per trajectory block");
print("<eigfrom> : lowest eigenvector index to take into account (e.g. 1)");
print("<eigto> : highest eigenvalue index to take into account (e.g. 10)");

#______________________________________________________________________________
## function psi: subspace overlap of eigensystem vector spaces A and B
##   formed by vectors with eigenvector index ranges 'irange' and 'jrange';
##   equation (2)
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

#______________________________________________________________________________
## function omega: covariance overlap;
##   equation (3)
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

#______________________________________________________________________________
## INPUT 
#______________________________________________________________________________
## PDB structure
pdb = read.pdb(args[1]);
print(pdb);

## DCD trajectory
dcd = read.dcd(args[2]);
print(dcd);

sBlock = as.numeric(ifelse(is.na(args[3]), 50, args[3]));
eigfrom = as.numeric(ifelse(is.na(args[4]), 1, args[4]));
eigto = as.numeric(ifelse(is.na(args[5]), 10, args[5]));

## eigenvalue range to consider; the same range for both blocks to compare
eigrange = eigfrom:eigto;
nEig = eigto - eigfrom + 1;
stopifnot(eigto <= dim(dcd)[2]);

#______________________________________________________________________________
## EIGENSYSTEMS of TRAJECTORY BLOCKS
#______________________________________________________________________________
## original size of trajectory
sTraj.orig = dim(dcd)[1];
## convert to integer multiplier
sTraj = floor(sTraj.orig / sBlock) * sBlock;
## create blocks of constant size 'sBlock'
block.startpos = seq(from = 1, to = sTraj, by = sBlock);
block.endpos = block.startpos + sBlock - 1;
block.pos = cbind(block.startpos, block.endpos);

## number of blocks in trajectory
nBlock = dim(block.pos)[1];

## lists of covariance matrices and eigensystems
covtraj = list(nBlock);
eigtraj = list(nBlock);

## compute covariance matrix and its eigensystem
for (i in 1:nBlock) {
	covtraj[[i]] = cov(dcd[(block.pos[i, 1]:block.pos[i, 2]), ]);
	eigtraj[[i]] = eigen(covtraj[[i]]);
}

#______________________________________________________________________________
## BLOCK PAIR OVERLAPS
#______________________________________________________________________________
## matrix to store overlaps between all block pairs
traj.overlap = matrix(0, nrow = nBlock, ncol = nBlock);

## matrix of omegaAB values of block pairs
for (i in 1:(nBlock-1)) {
	for (j in (i+1):nBlock) {
		traj.overlap[i, j] = omegaAB(eigtraj[[i]], eigrange, eigtraj[[j]], eigrange);
		traj.overlap[j, i] = traj.overlap[i, j];
	}
}

## Setting diagonal to '0' sets the colour range to that of the overlap values.
diag(traj.overlap) = 0;

#______________________________________________________________________________
## OUTPUT 
#______________________________________________________________________________
## show results as heatmap image
png("traj_overlap.png");
image(traj.overlap);
dev.off();

image(traj.overlap);


#===============================================================================

