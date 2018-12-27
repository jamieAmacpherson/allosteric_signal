#===============================================================================
# AlloHubMat
# Entropy metrics for a pair of character vectors
# Mutual Information, Finite Size Error, Joint Entropy, normalised MI
# See also our paper https://doi.org/10.1096/fj.11-190868.
# (C) 2018 Jens Kleinjung
#===============================================================================
## input data col1 and col2 should be character vectors
## might work with other data types, but that is untested
.mi = function(col1, col2) {
	## nput checks
	stopifnot(length(col1) > 1);
	stopifnot(is.character(col1));
	stopifnot(length(col2) > 1);
	stopifnot(is.character(col2));

	## marginal probabilities
	t1 = table(col1);
	p1 = t1 / sum(t1);
	t1.l = length(t1);

	t2 = table(col2);
	p2 = t2 / sum(t2);
	t2.l = length(t2);

	## vector of letter pair words
	pp.w = paste(c(col1), c(col2), sep = "");
	## vector of pairs ...
	t.pp = table(pp.w);
	t.pp.l = length(t.pp);
	##   ... and their marginal probabilities
	pp = t.pp / sum(t.pp);

	## vector of first letters of pairs ...
	pp.1 = substring(names(pp), 1, 1);
	##   ... and their marginal probabilities
	pp.p1 = p1[pp.1];

	## vector of second letters of pairs ...
	pp.2 = substring(names(pp), 2, 2);
	##   ... and their marginal probabilities
	pp.p2 = p2[pp.2];

	## Mutual Information
	MI = sum(pp * log(pp / (pp.p1 * pp.p2)));

	## Finite Size Error
	FSE = (t.pp.l - t1.l - t2.l + 1) / sum(length(col1), length(col2));

	## Joint Entropy
	JE = -sum(pp * log(pp));

	## normalised Mutual Information
	if (JE > 0) {
		nMI = (MI - FSE) / JE;
	} else {
		nMI = 0;
	}

	return(as.vector(c(MI, FSE, JE, nMI)));
}

#_______________________________________________________________________________
## returns random integer vector of length 'l' in value range [1, uplim]
.unif_int_range = function(l, uplim) {
  r_p = runif(l);
  r_i = as.integer((r_p * uplim) + 1);
  stopifnot(r_i > 0 & r_i < (uplim + 1));
  return(r_i);
}

#_______________________________________________________________________________
## computes MI, FSE, JE and nMI values between (all) column pairs
## careful: scales with square of sequence length
setMethod(f = "compute_mi_ali", signature = "matrix", definition = function(x) {
  #x@col_cmb = combn(length(x@refseq), 2);
  vars = as.numeric(names(x@pred_varpos_ver));
  stopifnot((length(vars) > 3) & (length(vars) < 30));
  nonvars = .unif_int_range(length(vars), length(x@refseq));
  x@varsnonvars = c(vars, nonvars);
  x@col_cmb = combn(x@varsnonvars, 2);
  cat(paste("MI for", dim(x@col_cmb)[2], "column pairs"));

  ## serial execution
  x@mi_mat = apply(x@col_cmb, 2, function(y) {
    .mi(x@ali[ , y[1]], x@ali[ , y[2]]);
  })

  ## parallel execution
  #cores = detectCores() - 1;
  #cl = makeCluster(cores);
  #x@mi_mat = parLapply(cl, x@col_cmb, 2, function(y) {
	#  .mi(x@ali[ , y[1]], x@ali[ , y[2]]);
	#})
  #stopCluster(cl);

  rownames(x@mi_mat) = c("MI", "FSE", "JE", "nMI");
  colnames(x@mi_mat) = apply(x@col_cmb, 2, function(x) paste(x, collapse = ":"));

  return(x);
})

#===============================================================================

