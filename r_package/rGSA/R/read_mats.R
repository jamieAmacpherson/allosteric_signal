#' Load the structural alphabet mutual information matrix trajectory.
#'
#' This function loads a series of symmetric matrices into the R environment.
#' It assumes that all mutual information matrices are symmetric and are sorted
#' according to their temporal order.
#'
#'
#' @param workingdir Path to the input file (eg. "~/pnurse/MD_dir/")
#' @param matrix.suffix Suffix for mutual information matrix used as a pattern identifier (eg. "*nMImat.out")
#' @return A list of mutual information matrices
#' @export



read_mats = function(workingdir, matrix.suffix) {

	print("READING MUTUAL INFORMATION MATRICES")

	details = file.info(list.files(pattern=matrix.suffix,
		full.names=TRUE,
		path = workingdir))

	details = details[with(details, order(as.POSIXct(mtime))), ]
	filenames = rownames(details)

	# read each of the mutual information matrices into the R environment
	datframe = lapply(filenames, read.table)

	# return a warning message if there are NA values in any of the matrices
	lapply(datframe, function(x) {
		if(is.na(x) == TRUE){
			writelines('WARNING: There is at least one NA value somewhere in your mutual information trajectory. \n All NA values will automatically assigned a zero value.')
		}
	}
	)

	# convert all NA values into 0
	datframe = lapply(datframe, function(x) { x[is.na(x)] <- 0; x})

	print("FINISHED READING MUTUAL INFORMATION MATRICES")

	# convert data frames into matrices
	listofmats = lapply(datframe, as.matrix)

	# return mutual information matrix trajectory as a list of matrices
	return(listofmats)

}

