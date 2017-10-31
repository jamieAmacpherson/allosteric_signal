#! /usr/bin/R

#===============================================================================
# comparison of POPSCOMP output for apo and fbp
#===============================================================================

## apo output 
inFileDirApo = "apo";
inFileNameApo = list.files(path = inFileDirApo, full.names = TRUE, pattern = 'dsasaMolecule$');

dataApo.l = lapply(inFileNameApo, function(x){read.table(file = x, header = TRUE)});
dataApo.df = do.call(rbind.data.frame, dataApo.l);
rownames(dataApo.df) = lapply(inFileNameApo, function(x) {substring(x, 5, 13)});
colnames(dataApo.df) = paste("APO", colnames(dataApo.df), sep = ".");

## fbp output 
inFileDirFbp = "fbp";
inFileNameFbp = list.files(path = inFileDirFbp, full.names = TRUE, pattern = 'dsasaMolecule$');

dataFbp.l = lapply(inFileNameFbp, function(x){read.table(file = x, header = TRUE)});
dataFbp.df = do.call(rbind.data.frame, dataFbp.l);
rownames(dataFbp.df) = lapply(inFileNameFbp, function(x) {substring(x, 5, 13)});
colnames(dataFbp.df) = paste("FBP", colnames(dataFbp.df), sep = ".");

## merge data frames
dataAll.df = merge(dataApo.df, dataFbp.df, by = "row.names");
write.table(dataAll.df, file = "popscomp_apo_fbp.out");

