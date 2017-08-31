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
# Start the clock!
ptm <- proc.time()
#______________________________________________________________________________


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
library("coop");

args = commandArgs(TRUE);
print("Usage: Rscript eigen_overlap.R <eigfrom> <eigto> <output prefix>"); 
print("<eigfrom> : lowest eigenvector index to take into account (e.g. 1)");
print("<eigto> : highest eigenvector index to take into account (e.g. 10)");
print("<output prefix> : output prefix");

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
		gpuMatMult(esA$vectors[ , irange[x]], t(esB$vectors[ , jrange[x]]));
	}));

	## using the 'abs' function to avoid negatives (from rounding errors)
	return (1 - sqrt(abs((ta - (2 * tb))) / ta));
}

#______________________________________________________________________________
## load Mutual information matrices 

readfiles = function() {
	details = file.info(list.files(pattern="*.out"))
	details = details[with(details, order(as.POSIXct(mtime))), ]
	filenames = rownames(details)
        datframe = lapply(filenames, read.table)
	return(lapply(datframe, as.matrix))
}
print("READING MUTUAL INFORMATION MATRICES")
nMImats = readfiles()
print("FINISHED READING MUTUAL INFORMATION MATRICES")

#______________________________________________________________________________
## EIGENSYSTEMS of mutual information matrices
#______________________________________________________________________________

## lowest eigenvector index of subspace
eigfrom = as.numeric(ifelse(is.na(args[1]), 1, args[1]));
## highest eigenvector index of subspace
eigto = as.numeric(ifelse(is.na(args[2]), 10, args[2]));

## eigenvalue range; the same range for both blocks to compare
eigrange = eigfrom:eigto;

nBlock = length(nMImats)

## lists of covariance matrices and eigensystems
covtraj = list(nBlock);
eigtraj = list(nBlock);

## compute covariance matrix and its eigensystem
print("COMPUTING COVARIANCE MATRICES")
covtraj = lapply(nMImats, covar);
print("COMPUTING EIGEN-SYSTEMS")
eigtraj = lapply(covtraj, eigen);


## matrix to store overlaps between all block pairs
traj.overlap = matrix(0, nrow = nBlock, ncol = nBlock);



#______________________________________________________________________________
## BLOCK PAIR OVERLAPS
#______________________________________________________________________________

## matrix of omegaAB values of block pairs
print("COMPUTING COVARIANCE OVERLAP OF MUTUAL INFORMATION MATRICE")
system.time(
for (i in 1:(nBlock-1)) {
	for (j in (i+1):nBlock) {
		traj.overlap[i, j] = omegaAB(eigtraj[[i]], eigrange, eigtraj[[j]], eigrange);
		traj.overlap[j, i] = traj.overlap[i, j];
	}
}
)


## CUDA version: matrix of omegaAB values of block pairs
#system.time(
#for (i in 1:(nBlock-1)) {
#	for (j in (i+1):nBlock) {
#		traj.overlap[i, j] = omegaAB_CUDA(eigtraj[[i]], eigrange, eigtraj[[j]], eigrange);
#		traj.overlap[j, i] = traj.overlap[i, j];
#	}
#}
#)

#______________________________________________________________________________
## OUTPUT 
#______________________________________________________________________________
## show results as heatmap image
diag(traj.overlap) = 0;

pdf(paste(args[3], "_cov_overlap.pdf", sep=""));
image.plot(traj.overlap);
dev.off();


#______________________________________________________________________________
## 2D KERNEL SMOOTHING
#______________________________________________________________________________
## set diagonal to '0' to limit overlap range to observed values
diag(traj.overlap) = 0;
## matrix values as vector
traj.overlap.v = as.vector(traj.overlap);
## matrix indices of vector values
nx = dim(traj.overlap)[1];
ny = dim(traj.overlap)[2];
x = rep(1:nx, ny);
y = rep(1:ny, each = nx);
grid.xy = as.data.frame(cbind(x,y));
## smooth vector values given grid indices, resulting in a list
traj.overlap.s = smooth.2d(traj.overlap.v, ind = grid.xy, nrow = nx, ncol = ny, theta = 2);

pdf(paste(args[3], "_cov_overlap_s.pdf", sep=""));
image.plot(traj.overlap.s);
dev.off();


#______________________________________________________________________________
## SPLIT OVERLAP MATRIX INTO SECTORS
#______________________________________________________________________________
## get the diagonal as a proxy for overlap values in the block 
sector.v = diag(traj.overlap.s$z);
pdf(paste(args[3], "_sectors.pdf", sep=""));
plot(sector.v, type = 's');
dev.off()

## split the diagonal into sectors, here done for each quantile
sector.liv = sapply(quantile(sector.v), function(x) x > sector.v);
## the 'FALSE' sectors at 75% should be useful
sector.sel.niv = which(! sector.liv[ ,"75%"]);
sector.idx = 1;
sector.idx.v = vector(length = length(sector.sel.niv));
sector.idx.v[1] = sector.idx
## count non-index-contiguous blocks
for (i in 2:length(sector.sel.niv)) {
	if (sector.sel.niv[i-1] != sector.sel.niv[i] - 1) {
		sector.idx = sector.idx + 1;
	}
	sector.idx.v[i] = sector.idx;
}

## translate block index vector into trajectory index vector
sector.sel.tiv = sector.sel.niv * sBlock;

## combine block index with trajectory index
sector.info = rbind(sector.idx.v, sector.sel.niv, sector.sel.tiv);
rownames(sector.info) = c("sector_ID", "block_idx", "traj_idx");
sector.info;
write.table(sector.info, file = paste(args[2], "_sectors.dat", sep = ""));

## save eigensystems and sector information for downstream trajectory comparisons
saveRDS(eigtraj, file = paste(args[2], "_eigtraj.RDS", sep = ""));
saveRDS(sector.info, file = paste(args[2], "_sectors.RDS", sep = ""));
saveRDS(eigrange, file = paste(args[2], "-eigrange.RDS", sep = ""));

#===============================================================================

# Stop the clock
proc.time() - ptm
