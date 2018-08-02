#! /usr/bin/python
##
## Perform Statistical Coupling Analysis on PKM2 multiple sequence
## alignment and compute the coupling of allosteric hubs identified
## from molecular dynamics simulations
##
##
##__________________________________________________________________
## Imports
import argparse
import os.path
import glob
import sys
from prody import *
from pylab import *
import matplotlib.pyplot as plt
from matplotlib import rcParams as pltparam

##__________________________________________________________________
## Plot styles
plt.style.use('seaborn-ticks')
pltparam.update({'font.size': 20})

##__________________________________________________________________
## Input parser
def is_valid_file(arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist! Use the --help flag for input options." % arg)
    else:
        return arg

parser = argparse.ArgumentParser(description='SCA analysis.')

parser.add_argument("-i", dest="fasta", required=True,
                    help="Multiple sequence alignment (FASTA)",
                    type=lambda x: is_valid_file(x))

args = parser.parse_args()

alignmentfile = args.fasta

##__________________________________________________________________
def sca_calc(MSAalignment, alignment_start, alignment_end):

	print '''
	Calculates the statistical coupling between all combinations of positions
	in a multiple sequence alignment.

	Inputs: 
	(1) Multiple sequence alignment.
	(2) First residue of the alignment.
	(3) last residue of the alignment.

	Outputs:
	(1) Heatmap plots of proximal statistical coupling scores between residues
	of the allosteric hubs, identified from molecular dynamics simulations.

	(2) Density histogram of the SCA matrix for the entire sequence alignment.

	(3) List of top-ranked SCA residues.
	'''

	# define the indices of alignment
	indices = list(range(alignment_start, alignment_end))

	# Parse the multiple sequence
	msa = parseMSA(MSAalignment)
	msa_refine = refineMSA(msa, label='PKM2_HUMAN', rowocc=0.8, seqid=0.98)

	# Calculate mutual information matrixmutinfo_norm = applyMutinfoNorm(mutinfo, entropy, norm='minent')
	entropy = calcShannonEntropy(msa_refine)
	mutinfo = buildMutinfoMatrix(msa_refine)

	# We can also apply normalization using applyMutinfoNorm() and correction using applyMutinfoCorr()
	# to the mutual information matrix based on references [Martin05] and [Dunn08], respectively
	mutinfo_norm = applyMutinfoNorm(mutinfo, entropy, norm='minent')
	mutinfo_corr = applyMutinfoCorr(mutinfo, corr='apc')

	plt.figure()
	showMutinfoMatrix(mutinfo_corr, clim=[0, 0.2], cmap='jet')
	plt.ylabel('Residue')
	plt.xlabel('Residue')
	plt.savefig('sca_pkm2.pdf', bbox_inches='tight')

	hubs = [121, 203, 244, 286, 304, 324, 355, 432, 489]
	hubs_dim = []

	for i in hubs:
		if i == 304:
			start = i - 0.5
			end = i + 4.5
			iterator = 5

		else:
			start = i - 0.5
			end = i + 3.5
			iterator = 4

		hubs_dim.append(range(i, i + iterator))

		plt.figure()
		showMutinfoMatrix(mutinfo_corr, clim=[0, 0.2], cmap='jet')
		plt.xlim(start, end)
		plt.ylim(start, end)
		plt.ylabel('Residue')
		plt.xlabel('Residue')
		plt.savefig('sca_%s_hub.pdf' %i, bbox_inches='tight')

	hubs_dim = np.asarray(hubs_dim)

	hubs_dim = np.hstack(hubs_dim)

	print(hubs_dim)

	hubs_mutinfo = mutinfo_corr[np.ix_(hubs_dim, hubs_dim)]
	print(np.shape(hubs_mutinfo))

	plt.figure()

	showMutinfoMatrix(hubs_mutinfo, clim=[0, 0.2], cmap='jet')
	plt.savefig('sca_pkm2_allhubs.pdf', bbox_inches='tight')

	## plot network of SCA couplings
	

	## plot histogram distribution of the SCA matrix
	
	# flatten SCA matrix into 1D array
	flat_sca = np.hstack(np.asarray(mutinfo_corr))
	mu = round(np.mean(flat_sca), 2)
	sigma = round(np.std(flat_sca), 2)

	plt.figure()
	n, bins, patches = plt.hist(flat_sca, 200, normed=1, facecolor='green', alpha=0.75)

	plt.xlabel('SCA score')
	plt.ylabel('Probability')
	plt.title(r'$\mu={0}, \sigma={1}$'.format(mu, sigma))
	plt.grid(True)
	plt.savefig('sca_pkm2_histogram.pdf', bbox_inches='tight')

	## Further analysis can also be done by rank ordering the matrix and analyzing
	# the pairs with highest mutual information or the most co-evolving residues.
	# This is done using calcRankorder(). A z-score normalization can also be applied
	# to select coevolving pairs based on a z score cutoff.


	print '''
	Rank-ordering performed by applying a y-score normalisation of all SCA couplings.

	Top 40 residues:
	'''
	rank_row, rank_col, zscore_sort = calcRankorder(mutinfo_corr, zscore=True)

	print asarray(indices)[rank_row[:20]], asarray(indices)[rank_col[:20]]
	

sca_calc(alignmentfile, 42, 395)
