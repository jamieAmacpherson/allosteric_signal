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

#______________________________________________________________________________
## load Mutual information matrices 

readMImat = function() {
	details = file.info(list.files(pattern="*nMImat.out", full.names=TRUE))
	details = details[with(details, order(as.POSIXct(mtime))), ]
	filenames = rownames(details)
        datframe = lapply(filenames, read.table)
	datframe = lapply(datframe, function(x) { x[is.na(x)] <- 0; x})
	return(lapply(datframe, as.matrix))
}

print("READING MUTUAL INFORMATION MATRICES")
nMImats = readMImat()
print("FINISHED READING MUTUAL INFORMATION MATRICES")


#______________________________________________________________________________
## load joint entropy matrices 

readjHmat = function() {
	details = file.info(list.files(pattern="*jHmat.out", full.names=TRUE))
	details = details[with(details, order(as.POSIXct(mtime))), ]
	filenames = rownames(details)
        datframe = lapply(filenames, read.table)
	datframe = lapply(datframe, function(x) { x[is.na(x)] <- 0; x})
	return(lapply(datframe, as.matrix))
}

print("READING ENTROPY MATRICES")
jHmats = readjHmat()
print("FINISHED READING MUTUAL INFORMATION MATRICES")

#______________________________________________________________________________
## calculate the total mutual information and joint entropy
## for a list of hubs 

hubs = c(24,70,178,197,205,238,288,357,362,390,418,440,492,512,311,306,307,310,356,413,489,491,515,62,305,237,287,417,63,347,377,490)

hubsa = hubs-13

MIout=list(hubsa)

for (i in hubsa){
	hubsa = lapply(nMImats, sum, [,i])
}
	

