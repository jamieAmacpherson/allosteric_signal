#! /usr/bin/R
#===============================================================================
# fit exponential
# Jens Kleinjung 2015
#===============================================================================

## command line arguments
args = commandArgs(TRUE);
stopifnot(! is.na(args[1]));
print(args);

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## RUN PARAMETERS
## 1. home path
PROJECTHOME = c("/home2/jkleinj/FRICTION");
# 2. domain name in upper case
DOM = args[1];
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

INPATH = paste(PROJECTHOME, "OUTPUT", "autocorrelation_runsNPT_ushort", DOM, sep = "/");
PROGRAMPATH = paste(PROJECTHOME, "PROGRAM", "fitexponential_runsNPT_ushort", sep = "/");
OUTPATH = paste(PROJECTHOME, "OUTPUT", "fitexponential_runsNPT_ushort", DOM, sep = "/");
DOMLC = tolower(DOM);

## prepare directories
#dir.create(PROGRAMPATH, showWarnings = TRUE, recursive = FALSE, mode = "0755");
dir.create(OUTPATH, showWarnings = TRUE, recursive = FALSE, mode = "0755");

## go to working directory
setwd(PROGRAMPATH);

#______________________________________________________________________________
## get the force files of all atoms
filenames = list.files(path = INPATH, full.names = FALSE, pattern = "^ac");

for (i in 1:length(filenames)) {
	print(i);
        INFILE = paste(INPATH, filenames[i], sep = "/");
        dat = scan(INFILE, quiet = TRUE);

        temp = data.frame(x = seq(length(dat)), y = dat);
	## catch fit errors with 'try' function
        dat.nls = try(nls(y ~ exp(- x / tau) + c, data = temp, start = list(tau = 10., c = 0.1)), silent = TRUE);

	## if fit successful, class of 'dat.nls' will be 'nls'
	##   otherwise 'try-error'
        if (class(dat.nls) == c("nls")) {

            ## write fitted 'tau' paramenter as result
            OUTFILE = paste(OUTPATH, "/", filenames[i], ".tau", sep = "");
            write.table(list(t(coef(dat.nls)), deviance(dat.nls)), file = OUTFILE, row.names = FALSE, col.names = FALSE);

            ## plot result
            OUTPLOT = paste(OUTPATH, "/", filenames[i], ".pdf", sep = "");
            pdf(OUTPLOT);
			par(mar=c(5, 5, 2, 2) + 0.1);
            plot(temp$x, temp$y, xlab = c("time / 100fs"), ylab = c("auto-correlation"), cex.lab = 2, cex.axis = 2, cex = 2);
            lines(temp$x, predict(dat.nls, list(x = temp$x)));
            dev.off();
        } else {
            cat("ERROR: exponential fit not converged");
        }
}

