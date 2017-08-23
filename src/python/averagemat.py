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
# Parse command line arguments 
#____________________________________________________________________________
def is_valid_file(arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist! Use the --help flag for input options." % arg)
    else:
	return arg

parser = argparse.ArgumentParser(description='Compute the average difference mutual information matrix.')


parser.add_argument('-c', type=float, nargs=1, 
                     help='Fixed value for the MP distribution', required=True)

parser.add_argument('-FBPpath', help= 'Path to FBP mutual information matrices (/path/to/directory)', type=str, action="store", required=True)

parser.add_argument('-APOpath', help= 'Path to Apo mutual information matrices (/path/to/directory)', type=str, action="store", required=True)

args = parser.parse_args()
#____________________________________________________________________________
# Compute the average MI matrices and take the difference matrix
#____________________________________________________________________________

def fbpmat(fbppath):

	path = fbppath + '/*.out'
	file_list = glob.glob(path)

	# sort files by date and time of creation
	file_list.sort(key=lambda x: os.path.getmtime(x))

	data = []
	for file_path in file_list:
		print "reading %s" %file_path
		data.append(
			np.loadtxt(file_path))
	
	# convert nan to zero values
	for i in data:
		where_are_NaNs = np.isnan(i)
		i[where_are_NaNs] = 0
	
	
	return(sum(data)/len(data))	


def apomat(apopath):

	path = apopath + '/*.out'
	file_list = glob.glob(path)

	# sort files by date and time of creation
	file_list.sort(key=lambda x: os.path.getmtime(x))

	data = []
	for file_path in file_list:
		print "reading %s" %file_path
		data.append(
			np.loadtxt(file_path))
	
	# convert nan to zero values
	for i in data:
		where_are_NaNs = np.isnan(i)
		i[where_are_NaNs] = 0
	
	return(sum(data)/len(data))	


def marcenkopasturpdf(x, c):
	# Marchenko Pastur density function for c > 1
	ub = (1 + math.sqrt(c))**2
	lb = (1 - math.sqrt(c))**2 
	mp = np.zeros(len(x))

	# Figure out indices where mp is to be calculated
	lbidx = np.where(x > lb)
	ubidx = np.where(x < ub)  
	a = lbidx[0][0]
	b = ubidx[-1][-1]
	xh = x[a:b+1]

	# MP distribution
	mp[a:b+1] = (((xh - lb)*(ub - xh))**0.5)/(2 * math.pi*c*xh)              
	return (lb, ub, mp)



def difmat(c):
	print 'Computing the Marcenko Pastur distribution function'
	
	# compute the average MI matrix for the holo structure
	fbpmatrix = fbpmat(args.FBPpath)
	np.savetxt('fbp_nMI.mat', fbpmatrix)

	# compute the average MI matrix for the apo structure
	apomatrix = apomat(args.APOpath)
	np.savetxt('apo_nMI.mat', apomatrix)

        # if the dimensions of the two matrices are different, resize the larger matrix
	# to fit the size of the smaller matrix. 
	if np.shape(fbpmatrix) != np.shape(apomatrix):
                print """##############################################################################
             WARNING: matrices are not of equal dimension
matrices will be resized to match the dimensions of the smallest input matrix.
This may lead to eraneous results. Proceed with caution.
##############################################################################"""
		dif = np.shape(fbpmatrix)[0]- np.shape(apomatrix)[0]

		if dif < 0:
			dif = dif * -1
			apomatrix = apomatrix[0:-dif, 0:-dif]

		else:
			fbpmatrix = fbpmatrix[0:-dif, 0:-dif]

		print "MP distribution with resized matrices"

        else:
                print "Computing MP distribution"

	matA = fbpmatrix - apomatrix
	np.savetxt('diff_mat.dat', matA)	

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
#	plt.plot(np.arange(0, 20, 0.01), mp, linewidth = 1, color='red')
	#plt.xlim(0,20)
	plt.ylim(-0.2,4) 
	plt.ylabel(r'P($\lambda$)')
	plt.xlabel(r'Eigenvalue, $\lambda$')
	plt.savefig('MMpdf.pdf', bbox_inches='tight')

difmat(args.c[0])
	
	


