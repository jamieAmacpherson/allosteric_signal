#! /usr/bin/R

#===============================================================================
# process all aligned pairs using POPSCOMP
#===============================================================================

library("bio3d");
## for atom selection mechanisms see:
## http://thegrantlab.org/bio3d/tutorials/structure-analysis

library("parallel");
## number of cores
nCore = detectCores() - 1;

## all PDB structures
inFileName = list.files(path = '.', full.names = FALSE, pattern = 'pdb$');

## output directory
outDir = "pops";

#_______________________________________________________________________________
## save A chain and split into domains
for (i in 1:length(inFileName)) {
	pdb = read.pdb(inFileName[i]);

	## we will use only chain A (the first of four protomers)
	chainA = atom.select(pdb, chain = "A");
	pdbChainA = trim.pdb(pdb, chainA);
	nameN = paste(substring(inFileName[i], 1, 6), "chainA", "pdb", sep = ".");
	write.pdb(pdbChainA, file = paste(outDir, nameN, sep = "/"));

	## N-terminal domain
	domN = atom.select(pdb, chain = "A", resno = 13:47);
	pdbN = trim.pdb(pdb, domN);
	nameN = paste(substring(inFileName[i], 1, 6), "domN", "pdb", sep = ".");
	write.pdb(pdbN, file = paste(outDir, nameN, sep = "/"));

	## A domain
	domA = atom.select(pdb, chain = "A", resno = c(48:116, 221:388));
	pdbA = trim.pdb(pdb, domA);
	nameA = paste(substring(inFileName[i], 1, 6), "domA", "pdb", sep = ".");
	write.pdb(pdbA, file = paste(outDir, nameA, sep = "/"));

	## B domain
	domB = atom.select(pdb, chain = "A", resno = 117:220);
	pdbB = trim.pdb(pdb, domB);
	nameB = paste(substring(inFileName[i], 1, 6), "domB", "pdb", sep = ".");
	write.pdb(pdbB, file = paste(outDir, nameB, sep = "/"));

	## C domain
	domC = atom.select(pdb, chain = "A", resno = 389:530);
	pdbC = trim.pdb(pdb, domC);
	nameC = paste(substring(inFileName[i], 1, 6), "domC", "pdb", sep = ".");
	write.pdb(pdbC, file = paste(outDir, nameC, sep = "/"));
}


#_______________________________________________________________________________
## POPS all structures
inFileName1 = as.list(list.files(path = 'pops', full.names = FALSE, pattern = 'pdb$'));

for (i in 1:length(inFileName1)) {
	#_______________________________________________________________________________
	## parallelisation
	## initiate cluster for parallel computation 
	clu = makeCluster(nCore);
	## make parallel functions see predefined variables
	clusterExport(clu, c("inFileName1"));

	#_______________________________________________________________________________
	## SAP alignment
	parLapply(clu, inFileName1, function(x) {
		system(paste("./pops --pdb pops/", x, sep = ""));
	});

	#_______________________________________________________________________________
	## release memory of parallelised structure
	stopCluster(clu);
};

#===============================================================================
