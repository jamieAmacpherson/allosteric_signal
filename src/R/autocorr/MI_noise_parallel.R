#! /usr/bin/R

#===============================================================================
# Compute noise statistics between successive matrices 
#===============================================================================

#______________________________________________________________________________
## LIBRARIES and FUNCTIONS
#______________________________________________________________________________
library("gtools");
library("parallel");

# number of cores
nCore <- detectCores() - 1;

#______________________________________________________________________________
## INPUT 
#______________________________________________________________________________
## MI matrices
inPath = c("./pkm2");
filenames = list.files(path = inPath, full.names = TRUE, pattern = 'lf_nMImat.out$');
filenames.sort = mixedsort(filenames);

tmp.df = read.table(filenames.sort[1]);
MI = rep(list(matrix(NA, nrow = dim(tmp.df)[1], ncol = dim(tmp.df)[2]),
			length(filenames.sort)));
for (i in 1:length(filenames.sort)) {
	print(paste(i, filenames.sort[[i]]));
	MI[[i]] = as.matrix(read.table(filenames.sort[i]));
}

#______________________________________________________________________________
## MI AUTOCORRELATION 
#______________________________________________________________________________
## autocorrelation and exponential fit
## MI.v: position-specific Mutual Information vector; i: iteration index
acMI = function(MI.v, i) {
	## autocorrelation
	MI_acf = acf(MI.v, plot = FALSE);

	## create data frame
	MI.df = data.frame(MI_acf$lag, MI_acf$acf);

	## exponential fit
	## catch fit errors with 'try' function
	MI.nls = try(nls(MI_acf.acf ~ exp(- MI_acf.lag / tau) + c, data = MI.df,
		start = list(tau = 1., c = 0.1)), silent = TRUE);

	## if fit successful, class of 'dat.nls' will be 'nls',
	##	otherwise it will be 'try-error'
	if (class(MI.nls) == c("nls")) {
		fitval[[i]] = c(coef(MI.nls), deviance(MI.nls));
	} else {
		fitval[[i]] = c(0, 0, 0);
	}
}

## upper triange lindices
rc = combn(1:dim(tmp.df)[1], 2);
## iteration indices
rc = rbind(rc, 1:dim(rc)[2]);
nRc = dim(rc)[2];

## data structure holding MI vectors of all positions
MI.vl = rep(list(vector(mode = "numeric", length = length(MI))), nRc);
for(i in 1:nRc) {
	MI.vl[[i]] = rapply(MI, function(x) x[rc[1,i], rc[2,i]]);
}

## data structure for fit values
fitval = rep(list(vector(mode = "numeric", length = 3)), nRc);
nIter = 0;

## initiate cluster for parallel computation 
clu = makeCluster(nCore);
## make parallel functions see predefined variables
clusterExport(clu, c("MI.vl", "nIter", "acMI", "fitval"));

invisible(parLapply(clu, MI.vl, function(x) { nIter <<- nIter + 1; acMI(x, nIter); }));

## run acMI function over all index values
## 'invisible' suppresses any output of the function evaluation
#invisible(sapply(1:nRc, function(x) { acMI(rc[1, x], rc[2, x], rc[3, x]); }));

saveRDS(fitval, "fitval.RDS");

#______________________________________________________________________________
stopCluster(clu);

#===============================================================================

