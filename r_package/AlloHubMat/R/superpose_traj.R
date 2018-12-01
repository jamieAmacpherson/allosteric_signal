#===============================================================================
# AlloHubMat
#' Superpose an MD trajectory to its starting coordinates, removing the rot-translational
#' degrees of freedom
#'
#' @param topol Topology file (pdb format)
#' @param traj Trajectory file (dcd format) 
#' @return superposed trajectory coordinates
#' @export
#===============================================================================

#_______________________________________________________________________________
## superpose trajectory frames
superpose_trj = function(topol, traj, traj.format = 'dcd'){

	if (identical(traj.format, "dcd")) {

		## select CA atoms for the trajectory frame superposition
		ca.inds <- bio3d::atom.select(topol, elety="CA");

		## fit the trajectory to the CA residues of the topology file
		xyz = bio3d::fit.xyz(fixed=topol$xyz, mobile=traj,
    	           fixed.inds=ca.inds$xyz,
    	           mobile.inds=ca.inds$xyz);

		## return the fitted trajectory coordinates
		return(xyz);

	} else{

		## dont support trajectory formats other than DCD for now
		stop(paste("trajectory extension", traj.format, 
			"not supported in this version of the package. Please convert your trajectory to DCD format. You can do this quickly and easily using the mdconvert tool (http://mdtraj.org/latest/mdconvert.html)"));
	}
} 

#===============================================================================

