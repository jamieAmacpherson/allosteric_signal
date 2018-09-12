#! /usr/bin/R

library("bio3d");

## read structure
pdb = read.pdb2("decaA.pdb");
## select CA atoms
ca.inds = atom.select(pdb, "calpha");
pdb$atom[ca.inds$atom, ];

## Example to test render_VMD.R
## artificial (random) data for 10 residues
cmb = combn(1:length(ca.inds$atom), 2);
dat = runif(dim(cmb)[2]);
cmb.dat = sapply(1:dim(cmb)[2], function(x) {
		cbind(cmb[1, x], pdb$atom[ca.inds$atom, ]$elety[cmb[1, x]],
			             pdb$atom[ca.inds$atom, ]$chain[cmb[1, x]],
			  cmb[2, x], pdb$atom[ca.inds$atom, ]$elety[cmb[2, x]],
			             pdb$atom[ca.inds$atom, ]$chain[cmb[1, x]],
			  dat[x]);
});

cmb.out = t(cmb.dat);
colnames(cmb.out) = c("resid1", "name1", "chain1", "resid2", "name2", "chain2", "corr");

## replace NAs in "chain"
cmb.out[is.na(cmb.out[ ,"chain1"]), "chain1"] = " ";
cmb.out[is.na(cmb.out[ ,"chain2"]), "chain2"] = " ";
cmb.out;

## select top 1/5 of correlations
sel.ix = cmb.out[ , "corr"] > 0.8;

## write correlations
write.table(cmb.out[sel.ix, ], "decaA_corrs.dat");
saveRDS(cmb.out[sel.ix, ], "decaA_corrs.RDS");

