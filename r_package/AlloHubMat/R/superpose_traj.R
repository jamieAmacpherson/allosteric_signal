#===============================================================================
# AlloHubMat
#' Superpose an MD trajectory to its starting coordinates,
#' removing the roto-translational degrees of freedom
#'
#' @param topol Topology file (pdb format)
#' @param traj Trajectory file (dcd format)
#' @return superposed trajectory coordinates
#' @export
#===============================================================================

#_______________________________________________________________________________
## superpose trajectory frames
superpose_trj = function(str, traj, traj.format = 'dcd') {
	if (identical(traj.format, "dcd")) {

		## select CA atoms for the trajectory frame superposition
		ca_inds = bio3d::atom.select(str, elety = "CA");

		## fit the trajectory to the CA residues of the topology file
		xyz = bio3d::fit.xyz(fixed = str$xyz, mobile = traj,
    	                   fixed.inds = ca_inds$xyz, mobile.inds = ca_inds$xyz);

		## return the fitted trajectory coordinates
		return(xyz);

	} else {
		## unsupported trajectory format
		stop(paste("trajectory extension", traj.format,
			"Unsupported trajectory format: Convert to 'dcd' using tools like mdconvert or VMD."));
	}
}

#===============================================================================
