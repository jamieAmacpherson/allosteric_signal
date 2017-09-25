#! /usr/bin/R

#===============================================================================
# align all pairs using SAP
#===============================================================================

library("parallel");
## number of cores
nCore = detectCores() - 1;

#_______________________________________________________________________________
## select structures
## 1. all PDB structures
#inFileName = list(1);

#inFileName[[1]] = list.files(path = '.', full.names = FALSE, pattern = 'pdb$');

## 2. by domain
#inFileName = list(4);

## domN domains
#inFileName[[1]] = list.files(path = '.', full.names = FALSE, pattern = 'domN.pdb$');
## domA domains
#inFileName[[2]] = list.files(path = '.', full.names = FALSE, pattern = 'domA.pdb$');
## domB domains
#inFileName[[3]] = list.files(path = '.', full.names = FALSE, pattern = 'domB.pdb$');
## domC domains
#inFileName[[4]] = list.files(path = '.', full.names = FALSE, pattern = 'domC.pdb$');

## 3. by domain and bound state
inFileName = list(8);

## domN APO domains
inFileName[[1]] = list.files(path = '.', full.names = FALSE,
	pattern = 'apo\\d{1}\\.domN\\.pdb$');
## domN FBP domains
inFileName[[2]] = list.files(path = '.', full.names = FALSE,
	pattern = 'fbp\\d{1}\\.domN\\.pdb$');

## domA APO domains
inFileName[[3]] = list.files(path = '.', full.names = FALSE,
	pattern = 'apo\\d{1}\\.domA\\.pdb$');
## domA FBP domains
inFileName[[4]] = list.files(path = '.', full.names = FALSE,
	pattern = 'fbp\\d{1}\\.domA\\.pdb$');

## domB APO domains
inFileName[[5]] = list.files(path = '.', full.names = FALSE,
	pattern = 'apo\\d{1}\\.domB\\.pdb$');
## domB FBP domains
inFileName[[6]] = list.files(path = '.', full.names = FALSE,
	pattern = 'fbp\\d{1}\\.domB\\.pdb$');

## domC APO domains
inFileName[[7]] = list.files(path = '.', full.names = FALSE,
	pattern = 'apo\\d{1}\\.domC\\.pdb$');
## domC FBP domains
inFileName[[8]] = list.files(path = '.', full.names = FALSE,
	pattern = 'fbp\\d{1}\\.domC\\.pdb$');

lapply(inFileName, function(i) {
	#_______________________________________________________________________________
	## structure pairs
	pair.cbn = combn(i, 2, simplify = FALSE);

	#_______________________________________________________________________________
	## parallelisation
	## initiate cluster for parallel computation 
	clu = makeCluster(nCore);
	## make parallel functions see predefined variables
	clusterExport(clu, c("pair.cbn"));

	#_______________________________________________________________________________
	## SAP alignment
	parLapply(clu, pair.cbn, function(x) { system(paste("sap", x[1], x[2])); });

	#_______________________________________________________________________________
	## release memory of parallelised structure
	stopCluster(clu);
});

#===============================================================================
