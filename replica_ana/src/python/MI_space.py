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
plt.style.use('seaborn-ticks')

#____________________________________________________________________________
# Parse the number of eigenmodes to be used in the calculation.
#____________________________________________________________________________
parser = argparse.ArgumentParser(description='Calculate the covariance overlap and cosine content of Mutual information matrices.')

parser.add_argument('nmodes', type=int, nargs=1, 
                     help='the number of eigenmodes used for the calculation')

parser.add_argument('--matrix', choices=['yes', 'no'], default='no', help='...')


args = parser.parse_args()

nmodes = args.nmodes[0]

#____________________________________________________________________________
# Read in mutual information matrix blocks generated by the GSAtools
# suite of programs.
#____________________________________________________________________________
path = '*.out'
file_list = glob.glob(path)
data = []
for file_path in file_list:
	print "reading mutual information matrix"
	data.append(
		np.loadtxt(file_path))



#____________________________________________________________________________
# Compute the cosine content:
#
# Mutual information matrices A and B are subjected to a spectral
# decomposition, to generate vectors of eigenvalues and matrices of 
# eigenvectors. The cosine content between the two matrices of eigenvectors
# is computed as a measure of matrix similarity.
#____________________________________________________________________________
def cosinecontent(matA, matB):
	
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
	
	return dab



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
	return omega


# compute the cosine content for the cartesian product of all MI matrices
def calcCC():
	veccosine = []
	linveccosine = []	

	for a in range(len(data)-1):
		b = a + 1
		print 'Calculating cosine content between blocks, using %s modes:' %nmodes
		print a, b
		linveccosine.append(cosinecontent(data[a], data[b]))

	np.savetxt('time_CosCont.dat', linveccosine)	

	# plot linear cosine content
	plt.plot(range(len(linveccosine)), linveccosine, color='black')	
    	plt.grid()
    	axes = plt.gca()
    	plt.ylabel(r'Cosine content, $\Psi_(A,B)$')
    	plt.xlabel(r'Simulation block')
    	plt.savefig('time_CosCont.pdf')
	
def calcCCmat():
	veccosine = []	
	for i, j in it.product(data, repeat=2):
		veccosine.append(cosinecontent(i, j))
	#np.savetxt('cosine_matrix.dat', veccosine)

	# reshape np array into matrix format
	cosinemat = np.reshape(veccosine, (len(data), len(data)))
	np.savetxt('cosine_matrix.dat', cosinemat)	

	plt.imshow(cosinemat, cmap='jet')
	plt.colorbar()
	#plt.clim(0,1)
	plt.savefig('cosineconent_mat.pdf')

plt.figure()	
calcCC()

# compute the covariance overlap for the cartesian product of all MI matrices
def calcCO():
	overlap = []
	
	for a in range(len(data)-1):
		b = a + 1
		print 'Calculating cosine content between blocks, using %s modes:' %nmodes
		print a, b
		overlap.append(covaroverlap(data[a], data[b]))
	np.savetxt('time_CovOverlap.dat', overlap)	
	
	# plot linear covariance overlap
	plt.plot(range(len(overlap)), overlap, color='red')	
    	plt.grid()
    	axes = plt.gca()
    	plt.ylabel(r'Covariance overlap, $\Omega_(A,B)$')
    	plt.xlabel(r'Simulation block')
    	plt.savefig('time_CovOver.pdf')

def calcCOmat():
	overlap = []
	for x, y in it.product(data, repeat=2):
		overlap.append(covaroverlap(x, y))
		
        # reshape np array into matrix format
        covarmat = np.reshape(overlap, (len(data), len(data)))

	np.savetxt('covar_overlap_matrix.dat', covarmat, delimiter='\t')
plt.figure()
calcCO()

if args.matrix == 'yes':
	plt.figure()
	calcCCmat()
else:
	exit()
