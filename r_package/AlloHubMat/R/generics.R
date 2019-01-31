#! /usr/bin/R

#===============================================================================
# Bayes variant detection
# (C) 2018 Jens Kleinjung
#===============================================================================

#_______________________________________________________________________________
#' sadata: An S4 class for Structural Alphabet data.
#' @slot fragment_letters: Letters forming the Structural Alphabet
#' @slot fragment_coordinates: Coordinates of prototype fragments
sadata <- setClass(
  "sadata",

  slots = c(
    fragment_letters = "vector",
    fragment_coordinates = "list",
    sa_trajectory = "matrix"
  )
)

#_______________________________________________________________________________
#' sector: An S4 class for trajectory sectors.
#' @slot contingency_pos_ver:  Contingency table of variant positions in vertical model
sector <- setClass(
  "sector",

  slots = c(
    contingency_pos_ver = "vector"
  )
)

#_______________________________________________________________________________
#' sector: An S4 class for hub residues.
#' @slot contingency_pos_ver:  Contingency table of variant positions in vertical model
hub <- setClass(
  "hub",

  slots = c(
    contingency_pos_ver = "vector"
  )
)

#_______________________________________________________________________________
## generic functions from read_str_traj.R
#_______________________________________________________________________________
#
#' Read input structure
#'
#' \code{read_str_file}
#'   reads a molecular structure file in PDB or GRO format.
#'
#' @param Input file name (including path) and format ("pdb" or "gro").
#' @return List of class "pdb".
#' @examples
#'   pdb_object = read_str_file(str_filename, str_format)
#'
setGeneric("read_str_file", function(x, ...) standardGeneric("read_str_file"));

#_______________________________________________________________________________
#
#' Read input trajectory
#'
#' \code{read_traj_file}
#'   reads a molecular trajectory file in DCD or XTC format.
#'
#' @param Input file name (including path) and format ("dcd" or "xtc").
#' @return .
#' @examples
#'   trajectory_object = read_traj_file(traj_filename, traj_format, start, end)
#'
setGeneric("read_traj_file", function(x, y, a, b, ...) standardGeneric("read_traj_file"));

#_______________________________________________________________________________
## generic functions from MI.R
#_______________________________________________________________________________
#
#' Compute Mutual Information matrix for pairs of alignment columns
#'
#' \code{compute_mi_ali}
#'   returns the metrics MI (Mutual Information), FSE (Finite Size Error), JE (Joint Entropy)
#'   and nMI (normalised Mutual Information) for column pairs of a given alignment.
#'   nMI is derived from the equation $nMI = (MI - FSE) / JE$.
#'
#' @param A character matrix.
#' @return A MI matrix.
#' @examples
#'   variBay_object = compute_mi_ali(matrix)
#'
setGeneric("compute_mi_ali", function(x, ...) standardGeneric("compute_mi_ali"));

#===============================================================================
