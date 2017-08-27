#! /usr/bin/R

#===============================================================================
# Compute the overlap of eigenspaces and their sub-spaces
# This CUDA version uses the gputools for accelerated matrix operations.
#===============================================================================
## Conformational Sampling and Dynamics of Membrane Proteins
##   From 10-Nanosecond Computer Simulations,
##   Faraldo-Gómez, Sansom et al.; DOI: 10.1002/prot.20257.
## see equations (2) and (3) therein
## It turns out that the normalised function 'omega' yields a more diverse range
##   of overlap values and therefore seems to be more suitable for splitting
##   trajectories into ergodic sectors (which appear as blocks with high
##   subspace overlap).

#______________________________________________________________________________
## LIBRARIES and FUNCTIONS
#______________________________________________________________________________

library("bio3d");
## Follow this installation instruction:
##   https://goatoftheplague.com/2016/12/08/installing-r-package-gputools-and-cuda-8-0-on-ubuntu-16-04/
## You might need to re-install the CUDA driver if the OS has been upgraded:
##    logout; CTRL-F1; su -; service lightdm stop; sh /usr/local/cuda*.run
library("gputools");
## http://www.image.ucar.edu/fields/
library("fields");

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
	#stopifnot(dim(esA$vectors[2]) >= range(irange)[2]);
	#stopifnot(dim(esB$vectors[2]) >= range(jrange)[2]);
	## assert index ranges are equally long
	#stopifnot(length(irange) == length(jrange));

	psiSum = sum(sapply(c(1:length(irange)), function(x) {
	  	esA$vectors[ , irange[x]] %*% esB$vectors[ , jrange[x]];
	}));
	return(psiSum / length(irange));
}

#______________________________________________________________________________
## function omega: covariance overlap;
##   equation (3)
omegaAB = function(esA, irange, esB, jrange) {
	## assert index range is within matrix dimension	
	#stopifnot(dim(esA$vectors[2]) >= range(irange)[2]);
	#stopifnot(dim(esB$vectors[2]) >= range(jrange)[2]);
	## assert index ranges are equally long
	#stopifnot(length(irange) == length(jrange));

	# sum terms of the equation are called here 'ta' and 'tb' */
	ta = sum(sapply(c(1:length(irange)), function(x) {
		esA$values[irange[x]] + esB$values[jrange[x]]; 
	}));

	tb = sum(sapply(c(1:length(irange)), function(x) {
		sqrt(esA$values[irange[x]] * esB$values[jrange[x]]) *
	  	esA$vectors[ , irange[x]] %*% esB$vectors[ , jrange[x]];		
	}));

	## using the 'abs' function to avoid negatives (from rounding errors)
	return (1 - sqrt(abs((ta - (2 * tb))) / ta));
}

omegaAB_CUDA = function(esA, irange, esB, jrange) {
	## assert index range is within matrix dimension	
	#stopifnot(dim(esA$vectors[2]) >= range(irange)[2]);
	#stopifnot(dim(esB$vectors[2]) >= range(jrange)[2]);
	## assert index ranges are equally long
	#stopifnot(length(irange) == length(jrange));

	# sum terms of the equation are called here 'ta' and 'tb' */
	ta = sum(sapply(c(1:length(irange)), function(x) {
		esA$values[irange[x]] + esB$values[jrange[x]]; 
	}));

	tb = sum(sapply(c(1:length(irange)), function(x) {
		sqrt(esA$values[irange[x]] * esB$values[jrange[x]]) *
		gpuMatMult(esA$vectors[ , irange[x]], esB$vectors[ , jrange[x]]);
	}));

	## using the 'abs' function to avoid negatives (from rounding errors)
	return (1 - sqrt(abs((ta - (2 * tb))) / ta));
}

#______________________________________________________________________________
## INPUT 
#______________________________________________________________________________
## PDB structure
args[1] = c("md1_ca.pdb"); 
pdb = read.pdb(args[1]);
print(pdb);

## DCD trajectory
args[2] = c("md1_rottrans_ca.dcd"); 
dcd = read.dcd(args[2]);
print(dcd);

## trajectory block size
sBlock = as.numeric(ifelse(is.na(args[3]), 50, args[3]));
## lowest eigenvector index of subspace
eigfrom = as.numeric(ifelse(is.na(args[4]), 1, args[4]));
## highest eigenvector index of subspace
eigto = as.numeric(ifelse(is.na(args[5]), 10, args[5]));
## eigenvalue range; the same range for both blocks to compare
eigrange = eigfrom:eigto;
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
system.time(
for (i in 1:(nBlock-1)) {
	for (j in (i+1):nBlock) {
		traj.overlap[i, j] = omegaAB(eigtraj[[i]], eigrange, eigtraj[[j]], eigrange);
		traj.overlap[j, i] = traj.overlap[i, j];
	}
}
)

## CUDA version: matrix of omegaAB values of block pairs
system.time(
for (i in 1:(nBlock-1)) {
	for (j in (i+1):nBlock) {
		traj.overlap[i, j] = omegaAB_CUDA(eigtraj[[i]], eigrange, eigtraj[[j]], eigrange);
		traj.overlap[j, i] = traj.overlap[i, j];
	}
}
)

#______________________________________________________________________________
## OUTPUT 
#______________________________________________________________________________
## show results as heatmap image
diag(traj.overlap) = 0;

png("traj_overlap.png");
image(traj.overlap);
dev.off();

image(traj.overlap);

#______________________________________________________________________________
## 2D KERNEL SMOOTHING
#______________________________________________________________________________
## matrix values as vector
diag(traj.overlap) = 0;

traj.overlap.v = as.vector(traj.overlap);
## matrix indices of vector values
nx = dim(traj.overlap)[1];
ny = dim(traj.overlap)[2];
x = rep(1:nx, ny);
y = rep(1:ny, each = nx);
grid.xy = as.data.frame(cbind(x,y));
## smooth vector values given grid indices, resulting in a list
traj.overlap.s = smooth.2d(traj.overlap.v, ind = grid.xy, nrow = nx, ncol = ny, theta = 2);

png("traj_overlap_s.png");
image(traj.overlap.s);
dev.off();

image(traj.overlap.s);

#______________________________________________________________________________
## SPLIT OVERLAP MATRIX INTO BLOCKS
#______________________________________________________________________________
## get the diagonal as a proxy for overlap values in the block 
bloc.v = diag(traj.overlap.s$z);
plot(bloc.v, type = 's');

## split the diagonal into sectors, here done for each quantile
bloc.liv = sapply(quantile(bloc.v), function(x) x > bloc.v);
## the 'FALSE' sectors at 75% should be useful
bloc.sel.niv = which(! bloc.liv[ ,"75%"]);
bloc.idx = 1;
bloc.idx.v = vector(length = length(bloc.sel.niv));
bloc.idx.v[1] = bloc.idx
## count non-index-contiguous blocks
for (i in 2:length(bloc.sel.niv)) {
	if (bloc.sel.niv[i-1] != bloc.sel.niv[i] - 1) {
		bloc.idx = bloc.idx + 1;
	}
	bloc.idx.v[i] = bloc.idx;
}

## combine block index with trajectory index
bloc.info = rbind(bloc.idx.v, bloc.sel.niv);
bloc.info;

#===============================================================================

