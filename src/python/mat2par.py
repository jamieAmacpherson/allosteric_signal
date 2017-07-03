#____________________________________________________________________________
# Compute convergence and similarity in mutual information space, between
# MD trajectories and over time within a single MD trajectory.
#
# Copyright 2017 Francis Crick Institute, King's College London
# and the Authors.
# 
# Author(s): Jamie A. Macpherson
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with GSAtools. If not, see <http://www.gnu.org/licenses/>.
#____________________________________________________________________________
# Imports
#____________________________________________________________________________
import argparse
import os.path
import glob
import sys
import numpy as np
import subprocess
import itertools as it
from numpy import linalg as LA
import math
import time
import progressbar
import matplotlib.pyplot as plt
from matplotlib import rcParams as pltparam
import sklearn.preprocessing as sk

plt.style.use('seaborn-ticks')
pltparam.update({'font.size': 20})

#____________________________________________________________________________
# Parse the number of eigenmodes to be used in the calculation.
#____________________________________________________________________________
def is_valid_file(arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist! Use the --help flag for input options." % arg)
    else:
	return arg

parser = argparse.ArgumentParser(description='Calculate the covariance overlap and cosine content of Mutual information matrices.')

parser.add_argument('nmodes', type=int, nargs=1, 
                     help='the number of eigenmodes used for the calculation')

parser.add_argument('-c', type=float, nargs=1, 
                     help='Fixed value for the MP distribution')

# the first argument is the trajectory file (.dcd) supplied after the -t flag
# the trajectory file is saved as an object with the variable args.dcdfile
parser.add_argument("-a", dest="mat1", required=True,
                    help="First Hermitian matrix",
                    type=lambda x: is_valid_file(x))

# the second argument is the topology file (.pdb) supplied after the -s flag
# this is saved an an obect with the variable args.pdbfile
parser.add_argument("-b", dest="mat2", required=True,
                    help="Second Hermitian matrix",
                    type=lambda x: is_valid_file(x))

args = parser.parse_args()

nmodes = args.nmodes[0]


#____________________________________________________________________________
# Compute the cosine content:
#
# Mutual information matrices A and B are subjected to a spectral
# decomposition, to generate vectors of eigenvalues and matrices of 
# eigenvectors. The cosine content between the two matrices of eigenvectors
# is computed as a measure of matrix similarity.
#____________________________________________________________________________
def cosinecontent(matA, matB):
	
	matA = np.loadtxt(matA)
	matB = np.loadtxt(matB)
	
	# kill if dimensions are different
	if np.shape(matA) != np.shape(matB):
		print "matrices must be of equal dimensions"
	else:
		print "Computing cosine content"
	
	# compute the eigenvalues and eigenvectors of matrix A
	eigvalA, eigvecA = LA.eig(matA) 
	
	# compute the eigenvalues and eigenvectors of matrix B
	eigvalB, eigvecB = LA.eig(matB)
	
	# initialize empty arrays to accept terms of cosine content measure
	cosinemult = []
	cosinesum = []
	
	for i in range(nmodes):
		# Vectors in an array such that the column v[:,i] is the eigenvector
		# corresponding to the eigenvalue w[i]
		cosinemult.append(eigvecA[:,i] * eigvecB[:,i])
		cosinesum.append(sum(cosinemult[i])**2)

	dab = sum(cosinesum) / nmodes
	
	print dab


#cosinecontent(args.mat1, args.mat2)

#____________________________________________________________________________
# Compute the covariance overlap:  
#
# Mutual information matrices A and B are subjected to a spectral
# decomposition, to generate vectors of eigenvalues and matrices of 
# eigenvectors. The covariance overlap is calculated by scaling the 
# cosine between eigenvectors A and B by the geometric mean of the respective
# eigenvalues.
#____________________________________________________________________________
def covaroverlap(matA, matB):
	matA = np.loadtxt(matA)
	matB = np.loadtxt(matB)
#
       # kill if dimensions are different 
        if np.shape(matA) != np.shape(matB):
                print "matrices must be of equal dimensions"
        else:
                print "Computing covariance overlap"
#
       # compute the eigenvalues and eigenvectors of matrix A
        eigvalA, eigvecA = LA.eig(matA)          
#
       # compute the eigenvalues and eigenvectors of matrix B
        eigvalB, eigvecB = LA.eig(matB)
#
       # initialize empty arrays to accept terms of covariance overlap measure
        valsum = []
        geomean = []
        geomeansum = []
        vecmult = []
        vecsum = []
	scaledcos=[]
#
 	#compute the sum of eigenvalues
	valsum = sum(eigvalA[0:nmodes] + eigvalB[0:nmodes])
	for i in range(nmodes):
#
      	# compute the sum of eigenvalues
      	# compute the geometric mean about two sets of eigenvalues
      		geomean.append(math.sqrt(eigvalA[i] * eigvalA[i]))
      	# determine the inner product of two eigenvector matrices
      		vecmult.append(eigvecA[:,i] * eigvecB[:,i])
      		vecsum.append(sum(vecmult[i])**2)
#
        # compute the vector scaled by the geometric mean of the 
        # eigenvalues	
	for x in range(nmodes):
		scaledcos.append(geomean[x] * vecsum[x])
#	
	scaledcossum = sum(scaledcos)
#
        omegain = valsum - round((2 * scaledcossum), 20)
        omega = 1 - (omegain/valsum)**0.5
#      
#
	print omega


#covaroverlap(args.mat1, args.mat2)


def marcenkopasturpdf(x, c):
	# Marchenko Pastur density function for c > 1
	ub = (1 + math.sqrt(c))**2
	lb = (1 - math.sqrt(c))**2 
	mp = np.zeros(len(x))
#
	# Figure out indices where mp is to be calculated
	lbidx = np.where(x > lb)
	ubidx = np.where(x < ub)  
	a = lbidx[0][0]
	b = ubidx[-1][-1]
	xh = x[a:b+1]
#
	# MP distribution
	mp[a:b+1] = (((xh - lb)*(ub - xh))**0.5)/(2 * math.pi*c*xh)              
	return (lb, ub, mp)


def plotmmpdf(mat, c):

	print 'Computing the Marcenko Pastur distribution function for %s' %mat
	matA = np.loadtxt(mat)
	
	nrow, ncol = np.shape(matA)
	N = nrow
	L = (N/c) 	

	# Scale matrix so that it has a zero mean and unit standard deviation
	#smatA = sk.scale(matA)	
	eigval,_ = LA.eig(matA) 
 	
	realeigval = np.array(eigval, dtype=float)

	weights = np.ones_like(realeigval)/float(len(realeigval))
	plt.hist(realeigval, bins=1000, normed=True, align='right', color='black')

	ln, un, mp = marcenkopasturpdf(np.arange(0, 20, 0.01), c)
	plt.plot(np.arange(0, 20, 0.01), mp, linewidth = 1, color='red')
	plt.xlim(0,20)
	plt.ylim(-0.2,4) 
	plt.ylabel(r'P($\lambda$)')
	plt.xlabel(r'Eigenvalue, $\lambda$')
	plt.savefig('MMpdf.pdf', bbox_inches='tight')



plotmmpdf(args.mat1, args.c[0])

