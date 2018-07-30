#' Read trajectory (XTC or DCD), return it as bio3d trajectory
#' Read structure (PDB or GRO), return it as bio3d structure
#'
#' @param workingdir Path to the input file (eg. "~/pnurse/MD_dir/")
#' @param name of the input trajectory and structure template file
#' @return bio3d structure trajectory
#' @export

library("bio3d");
library("streaMD");

#_______________________________________________________________________________
## read trajectory
## optionally XTC or DCD format, default is XTC
read_traj = function(workingdir, trajectory.name, first_frame, last_frame, traj.format = "xtc") {
	print("READING STRUCTURE TRAJECTORY");

	#workingdir = "/home/export/jens.kleinjung/jk.software/SAsuite/allosteric_signal/test_sys/deca_A_md/explicit/md";
	#trajectory.name = "decaA_1.xtc";

	## read XTC trajectory
	traj.path = paste(workingdir, trajectory.name, sep = "/");
	#first_frame = 1;
	#last_frame = 10;
	#traj.format = "xtc";
	if (identical(traj.format, "xtc")) {
		traj = loadxtc(traj.path, first_frame, last_frame);
		traj.bio3d = streaMD_to_bio3d(traj);
	} else if (identical(traj.format, "dcd")) {
		traj.bio3d = read.dcd(traj.path);
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
read_pdb = function(workingdir, structure.name, str.format = "pdb") {
	print("READING STRUCTURE");

	#workingdir = "/home/export/jens.kleinjung/jk.software/SAsuite/allosteric_signal/test_sys/deca_A_md/explicit/md";
	#structure.name = "decaA.pdb";
	str.format = "pdb";

	## read PDB structure
	str.path = paste(workingdir, structure.name, sep = "/");
	if (identical(str.format, "pdb")) {
		str.bio3d = read.pdb2(str.path);
	} else if (identical(str.format, "gro")) {
		str = loadgro(str.path);
		str.bio3d = gro2pdb(str);
	} else {
		stop(paste("structure extension", str.format, "not supported"));
	}

	print("FINISHED READING STRUCTURE");

	## return structure
	return(str.bio3d);
}

#===============================================================================


