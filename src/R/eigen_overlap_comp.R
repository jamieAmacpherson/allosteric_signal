#! /usr/bin/R

#===============================================================================
# Compare overlaps between sectors
#   Depends on RDS output from 'eigen_overlap_CUDA.R'
#===============================================================================

#______________________________________________________________________________
## LIBRARIES and FUNCTIONS
#______________________________________________________________________________

args = commandArgs(TRUE);
print("Usage: Rscript eigen_overlap_comp.R <system1> <system2>"); 
print("<system?> : molecular trajectories for which we have
	'<system?>-eigtraj.RDS', '<system?>-sectors.RDS' and '<system>-eigrange.RDS'");

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
		print(esA$values[irange[x]]);
		print(esA$values[jrange[x]]);
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
## INPUT 
#______________________________________________________________________________
## system 1
#args[1] = "md1_rottrans_ca.dcd";
system1 = args[1];
print(paste("System 1: ", system1));
eigtraj1 = readRDS(paste(args[1], "-eigtraj.RDS", sep = ""));
sectors1 = readRDS(paste(args[1], "-sectors.RDS", sep = ""));
eigrange1 = readRDS(paste(args[1], "-eigrange.RDS", sep = ""));

## system 2
#args[2] = "md1_rottrans_ca.dcd";
system2 = args[2];
print(paste("System 2: ", system2));
eigtraj2 = readRDS(paste(args[2], "-eigtraj.RDS", sep = ""));
sectors2 = readRDS(paste(args[2], "-sectors.RDS", sep = ""));
eigrange2 = readRDS(paste(args[2], "-eigrange.RDS", sep = ""));

stopifnot(eigrange1 == eigrange2);

#______________________________________________________________________________
## SECTOR PAIR OVERLAPS
#______________________________________________________________________________
## data structures for sector pair overlaps
nSector1 = max(sectors1["sector_ID", ]);
nSector2 = max(sectors2["sector_ID", ]);
traj.overlap.1_1 = matrix(0, nrow = nSector1, ncol = nSector1);
traj.overlap.2_2 = matrix(0, nrow = nSector2, ncol = nSector2);
traj.overlap.1_2 = matrix(0, nrow = nSector1, ncol = nSector2);
traj.overlap = list(traj.overlap.1_1, traj.overlap.2_2, traj.overlap.1_2);
names(traj.overlap) = c("sector_overlap.1_1", "sector_overlap.2_2", "sector_overlap.1_2");

## for all sector pairs
## 1_1
idx.comp = 1;
## all pairwise combinations of sector IDs
sectors1.l = split(sectors1["block_idx", ], sectors1["sector_ID", ]);
blocks1.unq = rapply(sectors1.l, function(x) head(x, 1));

## subspace overlap between two given (eig[idx]) eigensystems
for (i in 1:(length(blocks1.unq)-1)) {
	for (j in (i+1):length(blocks1.unq)) {
		traj.overlap[[idx.comp]][i, j] = omegaAB(eigtraj1[[blocks1.unq[i]]], eigrange1,
									 eigtraj1[[blocks1.unq[j]]], eigrange1);
		traj.overlap[[idx.comp]][j, i] = traj.overlap[[idx.comp]][i, j];
	}
}
diag(traj.overlap[[idx.comp]]) = 1;
png("traj_overlap.1_1.png");
image(traj.overlap[[idx.comp]]);
dev.off();

## 2_2
idx.comp = 2;
sectors2.l = split(sectors2["block_idx", ], sectors2["sector_ID", ]);
blocks2.unq = rapply(sectors2.l, function(x) head(x, 1));
for (i in 1:(length(blocks2.unq)-1)) {
	for (j in (i+1):length(blocks2.unq)) {
		traj.overlap[[idx.comp]][i, j] = omegaAB(eigtraj2[[blocks2.unq[i]]], eigrange2,
									 eigtraj2[[blocks2.unq[j]]], eigrange2);
		traj.overlap[[idx.comp]][j, i] = traj.overlap[[idx.comp]][i, j];
	}
}
diag(traj.overlap[[idx.comp]]) = 1;
png("traj_overlap.1_2.png");
image(traj.overlap[[idx.comp]]);
dev.off();

## 1_2
idx.comp = 3;
for (i in 1:(length(blocks1.unq)-1)) {
	for (j in (i+1):length(blocks2.unq)) {
		traj.overlap[[idx.comp]][i, j] = omegaAB(eigtraj1[[blocks1.unq[i]]], eigrange1,
									 eigtraj2[[blocks2.unq[j]]], eigrange2);
		traj.overlap[[idx.comp]][j, i] = traj.overlap[[idx.comp]][i, j];
	}
}
diag(traj.overlap[[idx.comp]]) = 1;
png("traj_overlap.1_2.png");
image(traj.overlap[[idx.comp]]);
dev.off();


#______________________________________________________________________________
## OUTPUT 
#______________________________________________________________________________
## save trajectory overlap object 
write.table(traj.overlap[[1]], file = "sector_eigen_overlap.1_1.dat");
write.table(traj.overlap[[2]], file = "sector_eigen_overlap.2_2.dat");
write.table(traj.overlap[[3]], file = "sector_eigen_overlap.1_2.dat");
saveRDS(traj.overlap, file = "sector_eigen_overlap.RDS");


#===============================================================================

