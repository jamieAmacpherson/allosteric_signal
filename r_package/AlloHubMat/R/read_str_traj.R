#===============================================================================
# AlloHubMat
#' Read trajectory (XTC or DCD), return it as bio3d trajectory
#' Read structure (PDB or GRO), return it as bio3d structure
#'
#' @param workingdir Path to the input file
#' @param name of the input trajectory and structure template file
#' @return bio3d structure trajectory
#' @export
#===============================================================================

#_______________________________________________________________________________
## read trajectory
## optionally XTC or DCD format, default is DCD
read_traj_file = function(workingdir, trajectory.name, traj.format = "dcd") {
	print("READING STRUCTURE TRAJECTORY");

	traj.path = paste(workingdir, trajectory.name, sep = "/");

	if (identical(traj.format, "dcd")) {
		traj.bio3d = bio3d::read.dcd(traj.path);
	} else if (identical(traj.format, "xtc")) {
		first_frame = 1;
		last_frame = 10;
		traj = streaMD::loadxtc(traj.path, first_frame, last_frame);
		traj.bio3d = streaMD::streaMD_to_bio3d(traj);
	} else {
		stop(paste("trajectory extension", traj.format, "not supported"));
	}

	print("FINISHED READING TRAJECTORY");

	## return trajectory
	return(traj.bio3d);
}

#_______________________________________________________________________________
## read structure
## optionally PDB or GRO format, default is PDB
read_pdb_file = function(workingdir, structure.name, str.format = "pdb") {
	print("READING STRUCTURE");

	## read PDB structure
	str.path = paste(workingdir, structure.name, sep = "/");
	if (identical(str.format, "pdb")) {
		str.bio3d = bio3d::read.pdb2(str.path);
	} else if (identical(str.format, "gro")) {
		str = streaMD::loadgro(str.path);
		str.bio3d = streaMD::gro2pdb(str);
	} else {
		stop(paste("structure extension", str.format, "not supported"));
	}

	print("FINISHED READING STRUCTURE");

	## return structure
	return(str.bio3d);
}

#===============================================================================

