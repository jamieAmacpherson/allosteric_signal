#' Plot the frequency distribution of the ergodic sector mutual information.
#' 
#' This function generates a frequency distribution plot of the mutual information contained by the ergodic sectors
#' identified by extract_sectors(), and calculates a significance cutoff to distinguish significant correlations
#' from the noise.
#'  
#' 
#'
#' @param sector.list A list of matrices, each is an averaged mutual information matrix for an ergodic sector
#' @param n.sigmas An integer number of standard deviations used in calculating the mutual information significance threshold (default value of 2)
#' @return A density histogram plot of the mutual information of the ergodic sectors
#' @export



sectors_MIfreq = function(sector.list, n.sigmas){

	# if the number of standard deviations is not specified by the user, assign default value of 2
	if(is.null(n.sigmas)){
		n.sigmas = 2
	}

	# if there are more than a single ergodic sector in the list, then average over the ergodic
	# sector mutual information matrices
	if(length(sector.list) > 1){
		
		# perform an element-wise averaging over the list of sector matrices
		# to determine the average mutual information matrix
		mat = Reduce("+", sector.list) / length(sector.list)

	} else {

		# if there is only a single ergodic sector identified in the MD trajectory, then use
		# the mutual information matrix of that single ergodic sector
		mat = sector.list[[1]]
	}


	# determine the average of the average mutual information matrix
	mu.mat = mean(mat)

	# determine the standard deviation of the average mutual information matrix
	sig.mat = sd(mat)

	# calculate the significance threshold
	threshold = mumat + (n.sigmas*sig.mat)

	# histogram of the average mutual information matrix
	mathist = hist(mat,
		plot=FALSE,
		breaks='Freedman-Diaconis')


	# plot the histogram
	plot(mathist$breaks[-1],
		mathist$counts,
		type='h',
		xlab = "nMI",
		ylab = expression(italic(f(nMI))),
		cex.lab = 2,
		cex.axis = 2,
		panel.first = grid())
	
	# line denoting the significance threshold 	
	abline(v=threshold,
		lty=2,
		col="red",
		lwd=2)
	
	# text
	text(threshold + 0.3, (max(mathist$density) / 2),
		paste("nMI > ", round(threshold, 3)),
		cex=2)
	
}