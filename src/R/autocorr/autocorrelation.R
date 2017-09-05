#! /usr/bin/R
#===============================================================================
# compute autocorrelation 
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

INPATH = paste(PROJECTHOME, "OUTPUT", "forcesplit_runsNPT_ushort", DOM, sep = "/");
PROGRAMPATH = paste(PROJECTHOME, "PROGRAM", "autocorrelation_runsNPT_ushort", sep = "/");
OUTPATH = paste(PROJECTHOME, "OUTPUT", "autocorrelation_runsNPT_ushort", DOM, sep = "/");
DOMLC = tolower(DOM);

## prepare directories
#dir.create(PROGRAMPATH, showWarnings = TRUE, recursive = FALSE, mode = "0755");
dir.create(OUTPATH, showWarnings = TRUE, recursive = FALSE, mode = "0755");

## go to working directory
setwd(PROGRAMPATH);

#______________________________________________________________________________
## get the force files of all atoms
filenames = list.files(path = INPATH, full.names = FALSE, pattern = "^Fsolvent");

for (i in 1:length(filenames)) {
	INFILE = paste(INPATH, filenames[i], sep = "/");
	dat = read.table(INFILE, header = FALSE);

	# 1. column is atom number
	atomnumber = dat[1, 1];

	# 2.-4.column are Fsolvent components x,y,z
	dat.x.acf = acf(as.vector(dat$V2), lag.max = 40);
	dat.y.acf = acf(as.vector(dat$V3), lag.max = 40);
	dat.z.acf = acf(as.vector(dat$V4), lag.max = 40);

	# vector length from components
	dat.x.v = as.vector(unlist(dat.x.acf$acf));
	dat.y.v = as.vector(unlist(dat.y.acf$acf));
	dat.z.v = as.vector(unlist(dat.z.acf$acf));

	dat.ac = sqrt(dat.x.v^2 + dat.y.v^2 + dat.z.v^2);
	dat.ac.norm = dat.ac / sqrt(3);

	# write result
	OUTFILE = paste(OUTPATH, "/ac.", atomnumber, ".dat", sep = "");
	write.table(dat.ac.norm, file = OUTFILE, row.names = FALSE, col.names = FALSE);
}

