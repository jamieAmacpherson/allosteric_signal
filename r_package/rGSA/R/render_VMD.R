#' render_VMD.R
#'
#' VMD 'Visualization State' coded in a TCL script.
#' 
#' Reads residue correlation list and writes a TCL script
#' that can be loaded into VMD via 'File'->'Load Visualization State...'.
#' Given the corresponding PDB file is present in the same directory,
#' the structure is rendered in 'NewCartoon' style with colour gradient
#' red->blue from N- to C-terminus. MI correlations between residues
#' are shown as lines colour-coded by the correlation strngth.
#'
#' @param corrfile pdbname
#'
#' @return
#'
#' @export

render_VMD = function(corrfile, pdbname) {
	## read correlation data
	corrs = readRDS(corrfile);

	## read structure
	#pdb = bio3d::read.pdb2(pdbname);
	## select CA atoms
	#ca.inds = atom.select(pdb, "calpha");

	## colour map
	cols = c("blue", "cyan", "orange", "red");

	## write VMD (TCL) script
	## open output file
	sink("show_corrs.vmd");

	## function to reset visualisation
	#cat("proc reset_viz {molid} {\n");
	#cat("## operate only on existing molecules\n");
	#cat("if {[lsearch [molinfo list] $molid] >= 0} {\n");
	#cat("    # delete all representations\n");
	#cat("    set numrep [molinfo $molid get numreps]\n");
	#cat("    for {set i 0} {$i < $numrep} {incr i} {\n");
	#cat("      mol delrep $i $molid\n");
	#cat("    }\n");
	#cat("    ## add new representations\n");
	#cat("    ## protein cartoon\n");
	#cat("    mol color Index\n");
	#cat("    mol representation NewCartoon\n");
	#cat("    mol selection all\n");
	#cat("    mol material Opaque\n");

	cat("## TCL script for VMD to show correlations as lines\n");

	## basic selections and initial view
	## load molecule
	cat("## load molecule\n");
	cat(paste("mol new", pdbname, "type pdb first 0 last -1\n"));

	## remove default representation
	cat("## remove default representation\n");
	cat("mol delrep 0 0\n");

	## selection of protein backbone
	cat("## selection of protein backbone\n");
	cat("set allostr [atomselect top \"protein and backbone\"]\n");

	## add new representation
	cat("## new protein representation\n");
	cat("mol color Structure\n");
	cat("mol representation Trace\n");
	cat("mol selection $allostr\n");
	cat("mol material Opaque\n");
    	cat("mol addrep 0\n");

	cat(sprintf("## show all %d correlations\n", dim(corrs)[1]));
	## show all correlations 
	for (i in 1:dim(corrs)[1]) {

		## determine colour, skip low correlations
		col.val = as.integer((as.numeric(corrs[i, "corr"]) + 0.1) * 4);
		if (col.val > 0) {
			cat(sprintf("## show correlation %s\n", i));

			## do not print empty chains
			if (corrs[1, "chain1"] == " ") {
				cat(sprintf("set corrsel0 [atomselect top \"resid %s and name CA\"]\n", corrs[i, "resid1"]));
			} else {
				cat(sprintf("set corrsel0 [atomselect top \"resid %s and name CA and chain %s\"]\n", corrs[i, "resid1"], corrs[i, "chain1"]));
			}
			if (corrs[1, "chain1"] == " ") {
				cat(sprintf("set corrsel1 [atomselect top \"resid %s and name CA\"]\n", corrs[i, "resid2"]));
			} else {
				cat(sprintf("set corrsel1 [atomselect top \"resid %s and name CA and chain %s\"]\n", corrs[i, "resid2"], corrs[i, "chain2"]));
			}
			cat("set pos0 [lindex [$corrsel0 get {x y z}] 0]\n");
			cat("set pos1 [lindex [$corrsel1 get {x y z}] 0]\n");
			cat(paste("draw color", cols[col.val], "\n"));
			cat("draw line $pos0 $pos1 width 2\n");
		}
	}

	## close function
	#cat("  }\n");
	#cat("}\n");

	## reset visualisation
	#cat("reset_viz 0");

	## close output file
	sink();
}

