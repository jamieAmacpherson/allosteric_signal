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
library("zoo");
library("denstrip");
library("matrixStats");
library("data.table");
library("MASS");

args = commandArgs(TRUE);
print("Usage: Rscript eigen_overlap.R <eigfrom> <eigto> <output prefix>"); 
print("<eigfrom> : lowest eigenvector index to take into account (e.g. 1)");
print("<eigto> : highest eigenvector index to take into account (e.g. 10)");
print("<block length> : time length of the trajectory, in ns (e.g. 400)");
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
	details = file.info(list.files(pattern="*nMImat.out", full.names=TRUE))
	details = details[with(details, order(as.POSIXct(mtime))), ]
	filenames = rownames(details)
        datframe = lapply(filenames, read.table)
	datframe = lapply(datframe, function(x) { x[is.na(x)] <- 0; x})
	return(lapply(datframe, as.matrix))
}

print("READING MUTUAL INFORMATION MATRICES")
nMImats = readfiles()
print("FINISHED READING MUTUAL INFORMATION MATRICES")

#______________________________________________________________________________
## load Mutual information matrices 

readentropyfiles = function() {
	details = file.info(list.files(pattern="*entropy.dat", full.names=TRUE))
	details = details[with(details, order(as.POSIXct(mtime))), ]
	filenames = rownames(details)
        datframe = lapply(filenames, read.table)
	datframe = lapply(datframe, function(x) { x[is.na(x)] <- 0; x})
	return(lapply(datframe, as.matrix))
}

#print("READING POSITIONAL ENTROPY DATA")
#entropydat = readentropyfiles()
#print("FINISHED READING POSITIONAL ENTROPY DATA")

#______________________________________________________________________________
## EIGENSYSTEMS of mutual information matrices
#______________________________________________________________________________

## lowest eigenvector index of subspace
eigfrom = as.numeric(ifelse(is.na(args[1]), 1, args[1]));
## highest eigenvector index of subspace
eigto = as.numeric(ifelse(is.na(args[2]), 10, args[2]));

## eigenvalue range; the same range for both blocks to compare
eigrange = eigfrom:eigto;

# number of trajectory blocks
nBlock = length(nMImats)

# trajectory block size (in ns)
sBlock = as.numeric(args[3]) / nBlock

## lists of covariance matrices and eigensystems
covtraj = list(nBlock);
eigtraj = list(nBlock);

## compute covariance matrix and its eigensystem
print("COMPUTING MATRIX INFORMATION CONTENT")
nMIsum = lapply(nMImats, sum);
nMIvar = lapply(nMImats, sd);
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
time = seq(from=0, to=args[3], length.out=nBlock)
diag(traj.overlap) = 0;

pdf(paste(args[4], "_cov_overlap.pdf", sep=""));
par(mar=c(5,5,2,2))
par(lab.cex=2)

image.plot(time, time, traj.overlap,
	legend.cex=2,
	cex=2,
	xlab="Time (ns)",
	ylab="Time (ns)",
	zlim=c(0, 0.7));

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

pdf(paste(args[4], "_cov_overlap_s.pdf", sep=""));
par(mar=c(5,5,2,2))

image.plot(time, time, traj.overlap.s$z,
        legend.cex=2,
        cex=2,
        xlab="Time (ns)",
        ylab="Time (ns)");

dev.off();

# write smoothed matrix to file
smoothedMI = replace(traj.overlap.s$z, traj.overlap.s$z == 0, NA)

write.table(as.matrix(smoothedMI),
		file=paste(args[4], "_cov_overlap_s.dat", sep=""),
		sep=" ",
		col.names=F,
		row.names=F)

#fit = try(fitdistr(as.matrix(smoothedMI), densfun="log-normal"))

#pdf(paste(args[4], "_overlap_dist.pdf", sep=""))
#par(mar=c(5,5,2,2))

#hist(as.matrix(smoothedMI),
#	breaks=20,
#	xlim=c(0,1),
#	xlab = "Covariance Overlap",
#	ylab = "Density",
#	prob=TRUE,
#	main="",
#	cex.axis=2,
#	cex.lab=2)

#lines(dnorm(as.matrix(smoothedMI),
#		fit$estimate[1],
#		fit$estimate[2]),
#		col="red")
#dev.off()

#fit

#______________________________________________________________________________
## SPLIT OVERLAP MATRIX INTO SECTORS
#______________________________________________________________________________
## get the diagonal as a proxy for overlap values in the block 
sector.v = diag(traj.overlap.s$z);
pdf(paste(args[4], "_sectors.pdf", sep=""));
par(mfrow=c(3,1))
plot(time, sector.v,
	type = 's',
	xlab = 'Time (ns)',
	ylab = expression(paste(Omega[ab])))
plot(time, unlist(nMIsum),
	type='s',
	xlab = 'Time (ns)',
	ylab = 'MI content')
plot(time, unlist(nMIvar),
	type='s',
	xlab = 'Time (ns)',
	ylab = 'MI variance')
dev.off()

# write the diagonal of the covariance overlap matrix to file
write.table(sector.v, file="cov_overlap_diag.dat", col.names=F, row.names=F)

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
write.table(sector.info, file = paste(args[4], "_sectors.dat", sep = ""));

## save eigensystems and sector information for downstream trajectory comparisons
saveRDS(eigtraj, file = paste(args[4], "_eigtraj.RDS", sep = ""));
saveRDS(sector.info, file = paste(args[4], "_sectors.RDS", sep = ""));
saveRDS(eigrange, file = paste(args[4], "_eigrange.RDS", sep = ""));
write.table(sector.info, file = paste(args[4], "_sectors.dat", sep=""));

#______________________________________________________________________________
## AVERAGE OVER CONTIGUOUS ERGODIC BLOCKS AND EXTRACT DISCRETE (AVERAGED) BLOCKS
#______________________________________________________________________________
extract.sectors = function(){
	sector.cont.ind = seqToIntervals(sector.info[2,])
	nSectors.cont = nrow(sector.cont.ind)
	sector.cont = list(nSectors.cont)

	dimx = nrow(nMImats[[1]])
	dimy = ncol(nMImats[[1]])

	sectors = list()
	sectorsentropy = list()

	## find contiguous sectors and add them to list "sectors"
	for (i in 1:nSectors.cont){
		from=sector.cont.ind[i,1]
		to=sector.cont.ind[i,2]
		print(from)
		print(to)
	
		## subset the contiguous sectors and find the element-
		## wise average of the sector
		# check whether we can average on the list directly
		submat = nMImats[from:to]
		Y = do.call(cbind, submat)
		Y = array(Y, dim=c(dim(submat[[1]]), length(submat)))
		sectors[[i]] = apply(Y, c(1, 2), mean, na.rm = TRUE)
		
		## write the element-averaged sectors to files	## put this outside the for loop
		write.table(sectors[[i]],
			file=paste("ergodic_sector_", i, sep=""),
			col.names=F, row.names=F,
			sep=" ")
		
#		# subset the entropy in the contiguous sectors
#		subentr = entropydat[from:to]
#		X = do.call(cbind, subentr)
#		X = array(X, dim=c(dim(subentr[[1]]), length(subentr)))
#		sectorsentropy[[i]] = apply(X, c(1, 2), mean, na.rm = TRUE)

		## write the element-averaged sector entopy to files  ## put this outside the for loop
#                write.table(sectorsentropy[[i]],
#                        file=paste("ergodic_sector_entropy_", i, sep=""),
#                        col.names=F, row.names=F,
#                        sep=" ")

	}

	## lists of covariance matrices and eigensystems
	sector.cov = list(nBlock);
	sector.eig = list(nBlock);

	## compute covariance matrix and its eigensystem
	print("COMPUTING COVARIANCE MATRICES OF ERGODIC SECTORS")
	sector.cov = lapply(sectors, covar);
	print("COMPUTING EIGEN-SYSTEMS OF ERGODIC SECTORS")
	sector.eig = lapply(sector.cov, eigen);
	
	## compute the covariance overlap between the ergodic sectors
	sector.overlap = matrix(0, nrow = nSectors.cont, ncol = nSectors.cont);
	
	print("COMPUTING COVARIANCE OVERLAP OF ERGODIC SECTORS")
	for (i in 1:(nSectors.cont-1)) {
		for (j in (i+1):nSectors.cont) {
			sector.overlap[i, j] = omegaAB(sector.eig[[i]], eigrange, sector.eig[[j]], eigrange);
			sector.overlap[j, i] = sector.overlap[i, j];
		}
	}

	## show results as heatmap image
	diag(sector.overlap) = 0;

	pdf(paste(args[4], "ergsector_cov_overlap.pdf", sep=""));
	par(mar = c(5,5,1,2))
	image.plot(sector.overlap,
		xlab = "Ergodic sector",
		ylab = "Ergodic sector",
		cex.lab = 2,
		cex.axis = 2);
	dev.off();
	
	
}
	  
extract.sectors()

#===============================================================================

# Stop the clock
proc.time() - ptm
