#! /usr/bin/R

#===============================================================================
# Bayes variant detection
# (C) 2018 Jens Kleinjung
#===============================================================================

#_______________________________________________________________________________
## S4 classes
#_______________________________________________________________________________
#' variBay: An S4 class for variant detection data.
#' @slot inpar:       parameters defining the alignment (variant) parameters
#' @slot runmode:     run parameters, e.g. cut-off values for classification
#' @slot refseq:      reference sequence
#' @slot ali:         alignment containing variants
#' @slot f_b_r:       positional frequency of reference base
#' @slot f_b_nr:      same of non-reference base with highest occurrence in column
#' @slot f_noise:     same of remaining bases plus gap
#' @slot post_ver:    posterior variant position probability from vertical model
#' @slot loglik_ver:  variant position log-likelihood from vertical model
#' @slot pred_varpos_ver: predicted variant position in vertical model
#' @slot varsnonvars: variant and non-variant positions for MI calculation
#' @slot col_cmb:     column pair combinations for MI calculation
#' @slot mi_mat:      MI matrix derived from column pairs
#' @slot post_hor:    posterior variant position probability from horizontal model
#' @slot loglik_hor:  variant position log-likelihood from horizontal model
#' @slot pred_pairpos_hor: predicted column-pair (MI) variant in horizontal model
#' @slot pred_varpos_hor: predicted variant position in horizontal model
#' @slot post_selpos: selection of variant positions: score and truth value
#' @slot post_selseq:  selection of variant sequences: score and truth value
variBay <- setClass(
  "variBay",

  slots = c(
  inpar = "list",
  runmode = "list",
  rnd_varpos = "numeric",
  rnd_varbas = "character",
  rnd_varseq = "numeric",
  refseq = "character",
  ali = "matrix",
  f_b = "list",
  f_b_r = "numeric",
  f_b_nr = "numeric",
  f_noise = "numeric",
  post_ver = "numeric",
  loglik_ver = "numeric",
  pred_varpos_ver = "numeric",
  varsnonvars = "numeric",
  col_cmb = "matrix",
  mi_mat = "matrix",
  class_hor = "logical",
  post_hor = "numeric",
  loglik_hor = "numeric",
  pred_pairpos_hor = "numeric",
  pred_varpos_hor = "numeric",
  post_selpos = "matrix",
  post_seq = "numeric",
  loglik_seq = "numeric",
  pred_seq = "numeric"
  )
)

#_______________________________________________________________________________
#' bema: An S4 class for benchmarking data.
#' @slot contingency_pos_ver:  Contingency table of variant positions in vertical model
#' @slot contingency_pos_hor:  Contingency table of variant positions in horizontal model
#' @slot contingency_pos_all:  Contingency table of variant positions in combined model
#' @slot contingency_seq_all:  Contingency table of variant sequences in combined model
bema <- setClass(
  "bema",

  slots = c(
    contingency_pos_ver = "vector",
    contingency_pos_hor = "vector",
    contingency_pos_all = "vector",
    contingency_seq_all = "vector"
  )
)

#_______________________________________________________________________________
## generic functions from MInorm.R
#_______________________________________________________________________________
#
#' Compute Mutual Information matrix for pairs of alignment columns (horizontal model)
#'
#' \code{compute_mi_ali}
#'   returns the metrics MI (Mutual Information), FSE (Finite Size Error), JE (Joint Entropy)
#'   and nMI (normalised Mutual Information) for column pairs of a given alignment.
#'   nMI is derived from the equation $nMI = (MI - FSE) / JE$.
#'
#' @param A variBay object.
#' @return A variBay object with filled-in slots 'col_cmb' (matrix) and 'mi_mat' (matrix).
#' @examples
#'   variBay_object = compute_mi_ali(variBay_object)
#'
setGeneric("compute_mi_ali", function(x, ...) standardGeneric("compute_mi_ali"));

#_______________________________________________________________________________
## generic functions from create_alignment.R
#_______________________________________________________________________________
#
#' Create sequence alignment
#'
#' \code{create_alignment}
#'   returns a random reference (RNA/DNA) sequence and
#'   sequence alignment given the input parameters.
#'   The random reference sequence of length given by 'inpar$lseq' is first stacked
#'   to a multiple alignment of identical sequences with a depth given by
#'   'inpar$nseq'. That alignment is then modified by introducing random SNVs
#'   with given positional frequency 'inpar$fvar_col' and variant sequence
#'   saturation 'inpar$far_row'. The final modification is a white noise
#'   sprinkler with a given per-base perturbation rate of 'inpar$fnoi'.
#'
#' @param A variBay object, 'inpar' parameters will be set on the fly.
#' @return A variBay object with filled-in slots 'inpar' (list),
#'   'refseq' (character vector) and 'ali' (character matrix).
#' @examples
#'   variBay_object = variBay();
#'   variBay_object = create_alignment(variBay_object)
#'
setGeneric("create_alignment", function(x, ...) standardGeneric("create_alignment"));

#_______________________________________________________________________________
## generic functions from compute_base_frequencies.R
#_______________________________________________________________________________
#
#' Compute base frequencies
#'
#' \code{compute_basefreqs}
#'   fills the slots 'f_b', 'f_b_r', 'f_b_nr' and 'f_noise'
#'   with base (plus gap) frequencies observed in the input alignment.
#'   The returned frequencies are vectors of observed frequencies [0,1],
#'   only 'f_b' is a list because of the potentially varying number of elements
#'   per position in shallow alignment.
#'
#' @param A variBay object.
#' @return A variBay object with filled-in slots 'f_b' (list), 'f_b_r' (vector),
#'   'f_b_nr' (vector) and 'f_noise' (vector).
#' @examples
#'   variBay_object = compute_basefreqs(variBay_object)
#'
setGeneric("compute_basefreqs", function(x, ...) standardGeneric("compute_basefreqs"));

#_______________________________________________________________________________
## generic functions from variBay.R
#_______________________________________________________________________________
#
#' Model variant likelihood per column position (vertical model)
#'
#' \code{model_vertical}
#'   uses pre-computed base frequencies (see function 'compute_basefreqs')
#'   and a signal-to-noise ratio to classify alignment columns as 'variant' or
#'   'non-variant'. Using that binary assignment, a Logistic Bayesian Regression
#'   is performed.
#'
#' @param A variBay object.
#' @return A fit object of class 'stanreg', 'glm' and 'lm' containing the
#'   regression results.
#' @examples
#'   fit_ver = model_vertical(variBay_object)

setGeneric("model_vertical", function(x, ...) standardGeneric("model_vertical"));

#_______________________________________________________________________________
#
#' Predict posterior probabilities per column position (vertical model)
#'
#' \code{predict_vertical}
#'   uses the posterior probabilities of the 'model_vertical'
#'   function to assertain variant positions. Generally positions with a
#'   probability >0.5 are considered variants.
#'
#' @param A variBay object.
#' @return A variBay object with filled-in slots 'post_ver' (vector),
#' 'loglik_ver' (vector) and 'pred_varpos_ver' (vector).
#' @examples
#'  variBay_object = predict_vertical(variBay_object)
#'
setGeneric("predict_vertical", function(x, ...) standardGeneric("predict_vertical"));

#_______________________________________________________________________________
#
#' Classify nMI vlaues as stemming from likely variant positions
#'
#' \code{make_class_mi}
#'   presuming a background MI distribution from non-variant positions,
#'   those positions with outlier-type high MI values are labelled as potential variants.
#'
#' @param A variBay object.
#' @return A variBay object with filled-in slots 'nmi_cutoff' (numeric) and
#'   'class_hor' (vector).
#' @examples
#'  variBay_object = make_class_mi(variBay_object)
#'
setGeneric("make_class_mi", function(x, ...) standardGeneric("make_class_mi"));

#_______________________________________________________________________________
#
#' Model variant likelihood derived from Mutual Information (horizontal model)
#'
#' \code{model_horizontal}
#'   uses pre-computed nMI values (see function 'compute_mi_ali')
#'   and a tentative variant classification (see function 'make_class_mi').
#'   Using that binary assignment, a Logistic Bayesian Regression is performed.
#'
#' @param A variBay object.
#' @return A fit object of class 'stanreg', 'glm' and 'lm' containing the
#'   regression results.
#'
#' @examples
#'  fit_hor = model_horizontal(variBay_object)
#'
setGeneric("model_horizontal", function(x, ...) standardGeneric("model_horizontal"));

#_______________________________________________________________________________
#
#' Predict posterior probabilities of Mutual Information pairs
#'
#' \code{predict_horizontal}
#'   uses the posterior probabilities of the 'model_horizontal'
#'   function to assertain variant positions. Generally positions with a
#'   probability >0.5 are considered variants.
#'
#' @param A variBay object.
#' @return A variBay object with filled-in slots 'post_hor' (vector),
#'   'loglik_hor' (vector) and 'pred_varpos_hor' (vector).
#' @examples
#'  variBay_object = predict_horizontal(variBay_object)
#'
setGeneric("predict_horizontal", function(x, ...) standardGeneric("predict_horizontal"));

#_______________________________________________________________________________
#
#' Combine results from vertical and horizontal models
#'
#' \code{combine_preds}
#'   combines results of the vertical and horizontal models to yield a
#'   comprehensive final result, i.e. predicted variant positions and
#'   variant sequences.
#'
#' @param A variBay object.
#' @return A variBay object with filled-in 'slots post_all' (vector),
#'   sel_varpos (vector) and sel_varseq (vector).
#' @examples
#'  variBay_object = combine_preds(variBay_object);
#'
setGeneric("combine_preds", function(x, ...) standardGeneric("combine_preds"));

#_______________________________________________________________________________
#
#' Model variant likelihood per sequence
#'
#' \code{model_sequence}
#'   uses the number (sum) of variants per sequence at predicted variant positions.
#'   The sums show a bimodal dsitribution and that is used to classify alignment
#'   sequences as 'variant' non-variant'. Using that binary assignment,
#'   a Logistic Bayesian Regression is performed.
#'
#' @param A variBay object.
#' @return A fit object of class 'stanreg', 'glm' and 'lm' containing the
#'   regression results.
#' @examples
#'   fit_seq = model_sequence(variBay_object)

setGeneric("model_sequence", function(x, ...) standardGeneric("model_sequence"));

#_______________________________________________________________________________
#
#' Predict posterior probabilities per sequence
#'
#' \code{predict_sequence}
#'   uses the posterior probabilities of the 'model_sequence'
#'   function to assertain variant positions. Generally positions with a
#'   probability >0.5 are considered variants.
#'
#' @param A variBay object.
#' @return A variBay object with filled-in slots 'post_seq' (vector),
#' 'loglik_seq' (vector) and 'pred_seq' (vector).
#' @examples
#'  variBay_object = predict_vertical(variBay_object)
#'
setGeneric("predict_sequence", function(x, ...) standardGeneric("predict_sequence"));

#_______________________________________________________________________________
## generic functions from variBay_bm.R
#_______________________________________________________________________________
#
#' Benchmark predicted variants
#'
#' \code{benchmark_preds}
#'   categorizes the predicted variants into the classes TP, FN, TN, FP
#'   if the true variant positions are given in the variBay slot 'rnd_varpos'.
#'
#' @param A variBay object and bema object.
#' @return A bema object.
#' @examples
#'   bema_object = benchmark_preds(variBay_object, bema_object)
#'
setGeneric("benchmark_preds", function(x, y, ...) standardGeneric("benchmark_preds"));

#_______________________________________________________________________________
## generic functions from read_alignment.R
#_______________________________________________________________________________
#
#' Benchmark predicted variants
#'
#' \code{read_alignment}
#'   reads input alignment in SAM format for given file path/name and
#'   assigns it as character matrix to passed variBay_object.
#'
#' @param A variBay object and input file path/name.
#' @return A bema object.
#' @examples
#'   bema_object = read_alignment(variBay_object, in_file_name)
#'
setGeneric("read_alignment", function(x, refseq, bam, bai, region, ...) standardGeneric("read_alignment"));

#===============================================================================
