#' Identify the allosteric hubs from the ergodic sector mutual information matrix  
#' 
#' This function computes the top n most signficant allosteric hub fragments. The function takes
#' the ergodic sector mutual information matrices from extract_sectors() and returns a list of the predicted top
#' most signficant allosteric hubs.
#'
#' The returned hubs are 'structural alphabet fragments' rather than amino acid residues. This level of abstraction
#' results from encoding the MD trajectory into strings of structural fragments. The per-residue mutual information
#' can be obtained by parsing the results of this function to resid_MIreconstr().
#'
#' 
#'  
#' 
#'
#' @param sector.list A list of matrices, each is an averaged mutual information matrix for an ergodic sector
#' @param n.hubs An integer number of the number of top most significant hubs.
#' @return A two-element list containing the top allosteric hubs and the associated mutual information couplings
#' between those allosteric hubs
#'
#' @export



identify_hubs = function(sector.list, n.hubs){


	# if the number of fragment hubs is not specified by the user, assign default value of 10
	if(is.null(n.hubs)){
		n.hubs = 10
	}

	if(n.hubs > ncol(sector.list[[1]])){
		writelines('WARNING \n Number of hubs selected exceeds total number of protein fragments.
			\n All of the protein fragments will be sorted according to their predicted
			contribution to the allosteric process. \n This may take a long time for large protein systems.')
	}


	# if there are more than a single ergodic sector in the list, then average over the ergodic
	# sector mutual information matrices
	if(length(sector.list) > 1){
		
		# perform an element-wise averaging over the list of sector matrices
		# to determine the average mutual information matrix
		mat = Reduce("+", sector.list) / length(sector.list)

		# do the same over the list of matrices to determine the element-wise
		# variance
		matvar = (apply(simplify2array(sector.list), 1:2, sd)) ** 2

		# generate a variance-weighted mutual information matrix, so that the couplings with
		# a low variance are up-weighted and those with a high-variance are down-weighted

		# transpose and translate the variances
		matvar.T = max(matvar) - matvar

		# variance-weight the mutual information matrix
		w.mat = mat * matvar.T

		# remove the lower triangle and diagonal elements of the 
		# mutual information matrix
		w.mat[lower.tri(w.mat)] = 0
		diag(w.mat) = 0

		# sort the mutual information according to the strength of the mutual information
		# couplings
		sorted.mat = sort(w.mat,
			decreasing = TRUE)

		# retrieve the matrix indices for the top n fragments
		# initialise a vector to contain the top n hubs
		top.hubs = c()

		# also capture the top mutual information couplings
		# initialise a data frame to capture the output
		top.couplings = data.frame()

		# loop through the desired number of top-ranked hubs and extract the coupled fragments which
		# are associated with the strongest coupling score (given by the highest mutual information)
		for(i in seq(from = 1, to = length(sorted.mat))){

			# determine the i, j indices of the sorted mutual information matrix
			tmp.hubs = which(mat == sorted.mat[i], arr.ind = TRUE)

			# append i and j fragments to dataframe <top.couplings>
			top.couplings[i, 1] = tmp.hubs[1]
			top.couplings[i, 2] = tmp.hubs[2]

			# append the indices to an intialised vector
			top.hubs = append(top.hubs, tmp.hubs)

		}

		# assign names to top.couplings dataframe
		names(top.couplings) = c('i', 'j')

		# extract a unique list of the top n number of hub fragments
		top.n.hubs = unique(top.hubs)[c(1:n.hubs)]

		# combine the unique list of hubs and the respective couplings 
		out.dat = list(top.n.hubs, top.couplings[c(1:n.hubs), ])

		# assign names to the list to be returned
		listnames = c('top.n.hubs', 'top.couplings')

		names(out.dat) = listnames

	}


	# if there is only a single ergodic sector identified in the MD trajectory
	else {

		# if there is only a single ergodic sector identified in the MD trajectory, then use
		# the mutual information matrix of that single ergodic sector
		mat = sector.list[[1]]


		# remove the lower triangle and diagonal elements of the 
		# mutual information matrix
		mat[lower.tri(mat)] = 0
		diag(mat) = 0


		# sort the mutual information according to the strength of the mutual information
		# couplings
		sorted.mat = sort(mat,
			decreasing = TRUE)

		# retrieve the matrix indices for the top n fragments
		# initialise a vector to contain the top n hubs
		top.hubs = c()

		# also capture the top mutual information couplings
		# initialise a data frame to capture the output
		top.couplings = data.frame()


		# loop through the desired number of top-ranked hubs and extract the coupled fragments which
		# are associated with the strongest coupling score (given by the highest mutual information)
		for(i in seq(from = 1, to = n.hubs)){


			# determine the i, j indices of the sorted mutual information matrix
			tmp.hubs = which(mat == sorted.mat[i], arr.ind = TRUE)

			# append i and j fragments to dataframe <top.couplings>
			top.couplings[i, 1] = tmp.hubs[1]
			top.couplings[i, 2] = tmp.hubs[2]

			# append the indices to an intialised vector
			top.hubs = append(top.hubs, tmp.hubs)

		}

		# assign names to top.couplings dataframe
		names(top.couplings) = c('i', 'j')

		# extract a unique list of the top n number of hub fragments
		top.n.hubs = unique(top.hubs)[c(1:n.hubs)]

		# combine the unique list of hubs and the respective couplings 
		out.dat = list(top.n.hubs, top.couplings[c(1:n.hubs), ])

		# assign names to the list to be returned
		listnames = c('top.n.hubs', 'top.couplings')

		names(out.dat) = listnames

	}

	# return a two-element list containing the top n hubs and the top couplings
	return(out.dat)

}


