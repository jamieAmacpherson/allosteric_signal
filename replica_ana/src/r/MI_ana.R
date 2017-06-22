#! /usr/bin/R

## load packages
library(fields)
library(TTR)
library(gplots)

# create list for all files ending with *.out
plot.matrices = function(){
	
	# Read all filenames ending with *.out into the R environment
	# these are all files containing the MI matrices at defined time points
	# during the simulation
	filenames = file.info(list.files(pattern=".out$", full.names=TRUE))
	filenames = filenames[with(filenames, order(as.POSIXct(mtime))), ]	
	
	# append the matrices into a list
	list.filenames=c()
	list.filenames = rownames(filenames)

	list.data <- list()
	for (x in 1:length(list.filenames)){
        	list.data[[x]] = as.matrix(read.table(list.filenames[x], header=FALSE))

		# Convert NA into zeros
		list.data[[x]][is.na(list.data[[x]])] <- 0.001
	}
	
	# add the names of the mutual information matrices to the list 
	names(list.data) <- list.filenames
	
	# generate sequences for i and j fragments
	i = seq(from=1, to=ncol(list.data[[1]]), by=1)
	j = seq(from=1, to=ncol(list.data[[1]]), by=1)


	# Plot the mutual information matrix as an image plot
#	png("image_plots.png")
#	
#	par(mfrow=c(2,2))
#	for (x in 1:length(list.data)) {
#	       	 title.name = gsub("\\D", "", list.filenames[[x]])
#		 image.plot(i, j, list.data[[x]],
#			main=paste(title.name, "ps", sep=" "))
#	}
#	
#	dev.off()
	
	# Make list.data a global variable so that the list can be accessed outside
	# of this function
	matrices <<- list.data
}

plot.matrices()


## Calculate the eigenvalues and eigenvectors of the mutual information
## matrices
eigen.matrix = function(){
	# Initialize a starting MI matrix with no intramolecular correlation
	start.mat = diag(nrow(matrices[[1]]))

	# Calculate eigenvalues and eigenvectors for the mutual information
	# matrices	
	for (x in 1:length(matrices)){
		matrices.eigen <<- eigen(matrices[[x]])
	}
}

## Select the number of modes to be included in cosine content and cosine
## overlap calculations
nmodes.cal = function(matrix.list){

}

## Function to calcualate the cosine of the angle between two eigenvectors
eigen.cosine = function(matrix.a, matrix.b){
	
	#Compute an eigenvalue decomposition of matrices A and B
	eigen.a = eigen(matrix.a)
	eigen.b = eigen(matrix.b)
	
	# Matrices A and B must be of the same dimensions
	if (ncol(matrix.a) != ncol(matrix.b)){
		stop("dimension mismatch")
	}

	# Compute the scaled dot product of two sets of eigenvectors from
	# matrices A and B.
	# The subspace overlap (d.ab) ranges from 0, when the eigenvector
	# subsets are orthogonal, to 1 when they are identical
	
	# Number of eigenmodes to be included in the calculation
	nmodes = 6 
        
	cosine.out = c()
	
	# loop through eigenmodes and determine the cosine content between eigenvector
	# matrices 
	cosine.out = t( t(eigen.a$vectors[,c(1:nmodes)]) * t(eigen.b$vectors[,c(1:nmodes)]))
	cosinesum = sum(colSums(cosine.out))
	d.ab = cosinesum/nmodes
#	for (i in 1:nmodes){
#		tmp = sum(eigen.a$vectors[,c(i)] * eigen.b$vectors[,c(i)])^2
#		cosine.out = append(cosine.out, tmp)
#	}
#
#	cosinesum = sum(cosine.out)
#	d.ab = cosinesum / nmodes
		
	return(d.ab)
}


## Calculate the covariance matrix overlap between two matrices
covariance.overlap = function(matrix.a, matrix.b){

	# Matrices A and B must be of the same dimensions
        if (ncol(matrix.a) != ncol(matrix.b)){
                stop("dimension mismatch")
        } 
	
	# Compute eigenvalue decomposition for matrices A and B
	eigen.a = eigen(matrix.a)
	eigen.b = eigen(matrix.b)	

	# Compute the upper term of equation (3) from Faraldo-Gomez, et al. Proteins
	# 2004 (Sampling in simulations of membrane proteins)
	
	# Number of eigenmodes to be included in the calculation
        nmodes = 6	

	# Sum over all eigenvalues of matrices and A and B
	omega.values = sum(eigen.a$values[c(1:nmodes)] + eigen.b$values[c(1:nmodes)])	
	
	upper.omega.vec = 2 * sum((sqrt(eigen.a$values[c(1:nmodes)] * eigen.b$values[c(1:nmodes)])) *
				(eigen.a$vectors[c(1:nmodes) ,c(1:nmodes)] %*%
					eigen.b$vectors[c(1:nmodes) ,c(1:nmodes)])^2)
 
	omega.ab = 1 - (((omega.values - upper.omega.vec)/
			omega.values)^0.5)

	return(omega.ab)

}

## Compute the cosine between eigenvector matrices of two mutual
## information matrices 
calc.eigen.cosine = function(matrix.list){
	
	tmp = c()
	
	for (x in 1:length(matrix.list)) {
		for (y in 1:length(matrix.list)) {
			
			tmp.out = eigen.cosine(matrix.list[[x]], matrix.list[[y]])
			
			tmp <- append(tmp, tmp.out)
			
		}
		
		calc.eigen.cosine.out = append(tmp, tmp.out)
	}	
	
	calc.eigen.cosine.out = head(calc.eigen.cosine.out, -1)
	
	cosine.matrix = matrix(calc.eigen.cosine.out,
				ncol = length(matrix.list),
				nrow = length(matrix.list))
	

	i = seq(from=1, to=length(matrix.list), by=1)
	j = seq(from=1, to=length(matrix.list), by=1)

	## Plot the mutual information matrix as an image plot
	pdf("cosine_matrix.pdf")
	
	image.plot(i, j, cosine.matrix,
		zlim = c(0,1),
		xlab = "Simulation block",
		ylab = "Simulation block",
		main = expression(paste(Psi[A:B])),
		cex.axis = 1.7, cex.lab = 1.7)
	
	
	dev.off()
	
	## Compute the cosine content for A_{t}:A_{t-1}.
	linear.eigen.cosine = c()
	
	# Loop over elements in list to compute covariance overlap between A_{t}:A_{t-1}
	for (i in 1:(length(matrix.list)-1)){
		t.1 = i
		t.2 = i + 1

		tmp.linear.cosine = (eigen.cosine(matrix.list[[t.1]], matrix.list[[t.2]]))
		
		linear.eigen.cosine <- append(linear.eigen.cosine, tmp.linear.cosine)
	}	
		
	# plot the covariance overlap between neighbouring blocks of the simulation 
	pdf("cosine_matrix_t.pdf")
	par(mar=c(5.1,5.1,5.1,2.1), mfrow=c(2,2))
	
	plot(linear.eigen.cosine,
		pch="",
		ylab=expression(paste(Psi[A[t]:A[t-1]])),
		xlab="t (blocks)",
		cex.lab=1.4, cex.axis=1.4,
		panel.first=grid())
	lines(linear.eigen.cosine)	



	# compute the first derivative of the cosine content to determine the
	# rate of change of covariance overlap (this should decay towards zero)
	roc.overlap = ROC(linear.eigen.cosine, type='discrete')

	# plot the first derivative of the cosine content
        plot(roc.overlap,
                pch="",
                ylab=expression(paste(Psi[A[t]:A[t-1]],"'")),
                xlab="t (blocks)",
                cex.lab=1.4, cex.axis=1.4,
                panel.first=grid())
	lines(roc.overlap)

	dev.off()	
			
}

calc.eigen.cosine(matrices)


## Compute the covariance overlap function for all pairs of blocks/sliding windows
## of the trajectory, using the mutual information matrix as an imput for the 
## analysis
calc.covariance.overlap = function(matrix.list){
	
	# Intialise a vector
	tmp = c()
	
	# loop over the list containing the mutual information matrices to select
	# for all possible combinations of matrices, including self-combinations
	for (x in 1:length(matrix.list)) {
		for (y in 1:length(matrix.list)) {
			
			# store the computed covariance overlap in a temporary vector
			tmp.out = covariance.overlap(matrix.list[[x]], matrix.list[[y]])
			
			# append the covariance overlap value into another temporary
			# vector
			tmp <- append(tmp, tmp.out)
			
		}
		
		# outside of the second nested loop, append the covariance overlap into
		# a final vector
		calc.eigen.cosine.out = append(tmp, tmp.out)
	}	

	# remove the final value from vector to generate a symmetric matrix vector	
	calc.eigen.cosine.out = head(calc.eigen.cosine.out, -1)
	
	# reshape the vector into a matrix
	cosine.matrix = matrix(calc.eigen.cosine.out,
				ncol = length(matrix.list),
				nrow = length(matrix.list))
	

	i = seq(from=1, to=length(matrix.list), by=1)
	j = seq(from=1, to=length(matrix.list), by=1)

	## Plot the mutual information matrix as an image plot
	pdf("covar_overlap_mat.pdf")
	
	image.plot(i, j, cosine.matrix,
		zlim = c(0,1),
		xlab = "Simulation block",
		ylab = "Simulation block",
		main = expression(paste("Covariance matrix overlap, ", Omega[A:B])),
		cex.axis = 1.7, cex.lab = 1.7)
	
	
	dev.off()

	## Compute the covariance overlap for A_{t}:A_{t-1}.
	linear.covar.overlap = c()
	
	# Loop over elements in list to compute covariance overlap between A_{t}:A_{t-1}
	for (i in 1:(length(matrix.list)-1)){
		t.1 = i
		t.2 = i + 1

		tmp.linear.covar = (covariance.overlap(matrix.list[[t.1]], matrix.list[[t.2]]))
		
		linear.covar.overlap <- append(linear.covar.overlap, tmp.linear.covar)
	}	
		
	# plot the covariance overlap between neighbouring blocks of the simulation 
	pdf("covar_overlap_t.pdf")
	par(mar=c(5.1,5.1,5.1,2.1), mfrow=c(2,2))
	
	plot(linear.covar.overlap,
		pch="",
		ylab=expression(paste(Omega[A[t]:A[t-1]])),
		xlab="t (blocks)",
		cex.lab=1.4, cex.axis=1.4,
		panel.first=grid())
	lines(linear.covar.overlap)	



	# compute the first derivative of the covariance overlap to determine the
	# rate of change of covariance overlap (this should decay towards zero)
	roc.overlap = ROC(linear.covar.overlap, type='discrete')
	
        plot(roc.overlap,
                pch="",
                ylab=expression(paste(Omega[A[t]:A[t-1]], "'")),
                xlab="t (blocks)",
                cex.lab=1.4, cex.axis=1.4,
                panel.first=grid())
	lines(roc.overlap)

	dev.off()	

	
	# Cluster the covariance overlap matrix to determine ergodic sectors within
	# the MD simulation
	
	pdf("ergodic_sectors.pdf")

	# Jet colour palette 
	jet.colours <-
		colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

	heatmap.2(cosine.matrix, trace="none",
			col=jet.colours,
			xlab = "Blocks",
			ylab = "Blocks")	
	
	dev.off()
}



#calc.covariance.overlap(matrices)



