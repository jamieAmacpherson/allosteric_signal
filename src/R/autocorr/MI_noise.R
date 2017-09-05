#! /usr/bin/R

#===============================================================================
# Compute noise statistics between successive matrices 
#===============================================================================

#______________________________________________________________________________
## LIBRARIES and FUNCTIONS
#______________________________________________________________________________
library("gtools");

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
nPos = dim(tmp.df)[1] * (dim(tmp.df)[1] - 1) / 2;
fitval = rep(list(vector(mode = "numeric", length = 5)), nPos);

## i: iteration index, r: row index, c: col index
acMI = function(i, r, c) {
	## autocorrelation given matrix position
	MI.v = rapply(MI, function(x) x[r, c]);
	MI_acf = acf(MI.v);

	## create data frame
	MI.df = data.frame(MI_acf$lag, MI_acf$acf);

	## exponential fit
	## catch fit errors with 'try' function
	MI.nls = try(nls(MI_acf.acf ~ exp(- MI_acf.lag / tau) + c, data = MI.df,
		start = list(tau = 1., c = 0.1)), silent = TRUE);

	## if fit successful, class of 'dat.nls' will be 'nls',
	##	otherwise it will be 'try-error'
	if (class(dat.nls) == c("nls")) {
		fitval[[i]] = c(r, c, coef(dat.nls), deviance(dat.nls));
	} else {
		fitval[[i]] = c(r, c, 0, 0, 0);
	}
}

#===============================================================================

