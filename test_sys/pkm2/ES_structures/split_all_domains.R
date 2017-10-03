#! /usr/bin/R

#===============================================================================
# split structures into separate domains 
#===============================================================================

library("bio3d");
## for atom selection mechanisms see:
## http://thegrantlab.org/bio3d/tutorials/structure-analysis

## all PDB pairs
filename = list.files(path = '.', full.names = FALSE, pattern = 'pdb$');

## output directory
outDir = "doms";

for (i in 1:length(filename)) {
	pdb = read.pdb(filename[i]);

	## we will use only chain A (the first of four protomers)
	## N-terminal domain
	domN = atom.select(pdb, resno = 13:47);
	pdbN = trim.pdb(pdb, domN);
	nameN = paste(filename[i], "domN", "pdb", sep = ".");
	write.pdb(pdbN, file = paste(outDir, nameN, sep = "/"));

	## A domain
	domA = atom.select(pdb, resno = c(48:116, 221:388));
	pdbA = trim.pdb(pdb, domA);
	nameA = paste(filename[i], "domA", "pdb", sep = ".");
	write.pdb(pdbA, file = paste(outDir, nameA, sep = "/"));

	## B domain
	domB = atom.select(pdb, resno = 117:220);
	pdbB = trim.pdb(pdb, domB);
	nameB = paste(filename[i], "domB", "pdb", sep = ".");
	write.pdb(pdbB, file = paste(outDir, nameB, sep = "/"));

	## C domain
	domC = atom.select(pdb, resno = 389:530);
	pdbC = trim.pdb(pdb, domC);
	nameC = paste(filename[i], "domC", "pdb", sep = ".");
	write.pdb(pdbC, file = paste(outDir, nameC, sep = "/"));
}

#===============================================================================
