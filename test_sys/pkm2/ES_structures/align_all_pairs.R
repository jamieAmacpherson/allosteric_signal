#! /usr/bin/R

#===============================================================================
# align all pairs using SAP
#===============================================================================

library("parallel");
## number of cores
nCore = detectCores() - 1;

## all PDB pairs
filenames = list.files(path = '.', full.names = FALSE, pattern = 'pdb$');
pair.cbn = combn(filenames, 2, simplify = FALSE);

## initiate cluster for parallel computation 
clu = makeCluster(nCore);
## make parallel functions see predefined variables
clusterExport(clu, c("pair.cbn"));

parLapply(clu, pair.cbn, function(x) { system(paste("sap", x[1], x[2])); });

## release memory
stopCluster(clu);

#===============================================================================
