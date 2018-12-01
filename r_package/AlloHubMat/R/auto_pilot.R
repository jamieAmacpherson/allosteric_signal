#==============================================================================
# AlloHubMat
#' Auto-pilot prediction of allosteric hub residues from a MD simulation.
#' 
#' This function reads the output of a molecular dynamics engine (in *.xtc format)
#' and a topology file (in *.pdb format), and runs through all the analysis functions
#' of AlloHubMat automatically. Where user-defined parameters are required,
#' default values are used.
#' Output files are stored in an automatically-generated directory.
#'
#' @param workingdir Path to the input file (eg. "~/sirpaulnurse/MD_dir/")
#' @param traj.length Length of the molecular dynamics simulation (in ns)
#' @return Predicted allosteric hub residues
#' @export 
#==============================================================================

auto_pilot = function(workingdir, traj.length){

	# Load the structural alphabet mutual information matrix trajectory
	mi.mats = read_mats(workingdir, matrix.suffix = '*.out')

	# Compute the eigensystem for the list of mutual information matrices
	eigensys = comp_eigensystem(mi.mats)

	# Compute the covariance overlap between each of the trajectory blocks from an input list
	# of eigensystems.
	overlap.mat = block_overlap(eigensys)

	# Smooth the matrix of covariance overlaps using a kernal smoothing algorithm.
	overlap.mat.s = matrix_smooth(overlap.mat)

	# Automated detection of ergodic sectors from a list of time-averaged mutual information matrices.
	erg.sec.info = detect_sectors(overlap.mat.s, traj.length)

	# Extract the ergodic sectors from a trajectory of mutual information matrices.
	erg.secs = extract_sectors(erg.sec.info, mi.mats)

	# Identify the allosteric hubs from the ergodic sector mutual information matrix  
	allosteric.hubs = identify_hubs(erg.secs) 

	## plotting functions
	
	# 2D ergodic sector plot
	pdf('2d_ergsector_plt.pdf')
	sectors_2dplot(overlap.mat.s, traj.length)
	dev.off()

	writeLines('2D plot generated identifying ergodic sectors within the MD simulation.')

	system('sleep -n 5')

	# 3D ergodic sector plot
	pdf('3d_ergsector_plt.pdf')
	sectors_3dplot(overlap.mat.s, traj.length)
	dev.off()

	writeLines('3D plot generated identifying ergodic sectors within the MD simulation.')

	system('sleep -n 5')

	return(allosteric.hubs$top.n.hubs)

}

#==============================================================================

