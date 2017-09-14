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
filenames = list.files(path = inPath, full.names = TRUE,
	pattern = '^[0-9][0-9][0-9][0-9][0-9].lf_nMImat.out$');
filenames.sort = mixedsort(filenames);

tmp.df = read.table(filenames.sort[1]);
MI = rep(list(matrix(NA, nrow = dim(tmp.df)[1], ncol = dim(tmp.df)[2]),
			length(filenames.sort)));
#for (i in 1:length(filenames.sort)) {
for (i in 1:length(filenames.sort)) {
	print(paste(i, filenames.sort[[i]]));
	MI[[i]] = as.matrix(read.table(filenames.sort[i]));
}

#______________________________________________________________________________
## MI AUTOCORRELATION 
#______________________________________________________________________________
## autocorrelation and exponential fit
## MI.v: position-specific Mutual Information vector; i: iteration index
acMI = function(MI.v) {
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
		c(coef(MI.nls), deviance(MI.nls));
	} else {
		c(0, 0, 0);
	}
}

#______________________________________________________________________________
## upper triange lindices
#rc = combn(1:dim(tmp.df)[1], 2);
## up to 500 works fine
rc = combn(1:500, 2);
## iteration indices
rc = rbind(rc, 1:dim(rc)[2]);
nRc = dim(rc)[2];

## data structure holding MI vectors of all positions
## that is following the MI signal of each position over time
MI.vl = rep(list(vector(mode = "numeric", length = length(MI))), nRc);
for(i in 1:nRc) {
	MI.vl[[i]] = rapply(MI, function(x) x[rc[1,i], rc[2,i]]);
}

## data structure for fit values
fitval = rep(list(vector(mode = "numeric", length = 3)), nRc);

#______________________________________________________________________________
## initiate cluster for parallel computation 
clu = makeCluster(nCore);
## make parallel functions see predefined variables
clusterExport(clu, c("MI.vl", "acMI", "fitval"));

## perform all fits
fitval = parLapply(clu, MI.vl, function(x) { acMI(x); });
## save results
saveRDS(fitval, "fitval.RDS");
## release memory
stopCluster(clu);

#______________________________________________________________________________
## evaluate fit results
## remove zero elements
fitval.nozero = fitval[lapply(fitval, sum) > 0];
## extract taus
tau = rapply(fitval.nozero, function(x) x[1]);
## tau distribution
hist(tau, breaks = 200, xlim = c(0, 5));


#===============================================================================

