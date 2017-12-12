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

parser.add_argument('-c', type=float, nargs=1, 
                     help='Fixed value for the MP distribution')
args = parser.parse_args()




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
	# convert nan to zero values
	for i in matA:
		where_are_NaNs = np.isnan(i)
		i[where_are_NaNs] = 0

	X = np.sqrt(0.5)*(np.randn(N, L) + 1j*np.randn(N, L))
	Wc = np.dot(X, conj(X.T))/L
	D, U = LA.eig(Wc)	
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

