#===============================================================================
# AlloHubMat
#' Read structure (PDB) and return it as bio3d structure
#' Read trajectory (DCD) and return it as bio3d trajectory
#'
#' @param workingdir Path to the input file
#' @param name of the input trajectory and structure template file
#' @return bio3d structure trajectory
#' @export
#===============================================================================

#_______________________________________________________________________________
## read structure
.get_str = function(input_str_file, str_format = "pdb") {
	print("Reading input structure\n");

	## PDB or error
	if (identical(str_format, "pdb")) {
		str_bio3d = bio3d::read.pdb2(input_str_file);
	} else {
		stop(paste("structure extension", str_format, "not supported"));
	}

	## return structure
	return(str_bio3d);
}

#_______________________________________________________________________________
setMethod(f = "read_str_file", signature = c("character"), definition = function(x, y) {
  .get_str(x, y);
})

#_______________________________________________________________________________
#_______________________________________________________________________________
## read trajectory
.get_traj = function(input_traj_file, traj_format = "dcd", first_frame = 1, last_frame = 10) {
	print("Reading input structure");

  ## DCD or error
	if (identical(traj_format, "dcd")) {
		traj_bio3d = bio3d::read.dcd(input_traj_file);
	} else {
		stop(paste("trajectory extension", traj_format, "not supported"));
	}

	## return trajectory
	return(traj_bio3d);
}

#_______________________________________________________________________________
setMethod(f = "read_traj_file", signature = c("character", "character", "numeric", "numeric"), definition = function(x, y, a, b) {
  .get_traj(x, y, a, b);
})

#===============================================================================
