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
from prody import *

plt.style.use('seaborn-ticks')
pltparam.update({'font.size': 20})


#____________________________________________________________________________
# Parse command line arguments____________________________________________________________________________
parser = argparse.ArgumentParser(description='Calculate the covariance overlap and cosine content of Mutual information matrices.')

parser.add_argument('nmodes', type=int, nargs=1, 
                     help='the number of eigenmodes used for the calculation.')

parser.add_argument('--COsum', choices=['yes', 'no'], default='no', help='(no/yes) Calculate the cumulative sum of the covariance overlap.')

parser.add_argument('--CO', choices=['yes', 'no'], default='no', help='(no/yes) Calculate the covariance overlap.')

parser.add_argument('--CCmat', choices=['yes', 'no'], default='no', help='(no/yes) Calculate the cosine content matrix.')

parser.add_argument('--COmat', choices=['yes', 'no'], default='no', help='(no/yes) Calculate the covariance overlap matrix.')

parser.add_argument('--PDBload', choices=['yes', 'no'], default='no', help='(no/yes) Parse PDB files.')

args = parser.parse_args()

nmodes = args.nmodes[0]



#____________________________________________________________________________
# Read in mutual information matrix blocks generated by the GSAtools
# suite of programs.
#____________________________________________________________________________
# Read in mutual information matrices
path = '*.out'
file_list = glob.glob(path)
data = []
for file_path in file_list:
	print "reading mutual information matrix"
	data.append(
		np.loadtxt(file_path))

# convert nan to zero values
for i in data:
	where_are_NaNs = np.isnan(i)
	i[where_are_NaNs] = 0

del i

# Read in PDB coordinates corresponding to the mutual information matrices
if args.PDBload == 'yes':
	pathpdb = '*.pdb'
	file_list = glob.glob(pathpdb)
	pdbdata = []
	for pdbfile_path in file_list:
		print "reading PDB coordinates"
		# read coordinates
		pdb = parsePDB(pdbfile_path)

		# remove the last three C $\alpha$ atoms
		pdb = pdb[0:(len(pdb)-3)].getCoords()
		pdbdata.append(pdb)


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
       # include the desired number of eigenmodes for matrix A
	eigvalA = eigvalA[0:nmodes]
	eigvecA = eigvecA[:,0:nmodes]
#
       # compute the eigenvalues and eigenvectors of matrix B
        eigvalB, eigvecB = LA.eig(matB)
#
       # include the desired number of eigenmodes for matrix B
	eigvalB = eigvalB[0:nmodes]
	eigvecB = eigvecB[:,0:nmodes]
#
#
	dotAB = np.dot(eigvecA.T, eigvecB)**2
#
	outerAB = np.outer(eigvalA**0.5, eigvalB**0.5)
#	
	diff = (np.sum(eigvalA.sum() + eigvalB.sum()) - 2 * np.sum(outerAB * dotAB))
#
	if diff < 0:
		diff = 0
	else:
		diff = diff**0.5 
	omega = 1 - diff / np.sqrt(eigvalA.sum() + eigvalB.sum())
#
	return omega










#def covaroverlap(matA, matB):
#
       # kill if dimensions are different 
#        if np.shape(matA) != np.shape(matB):
#                print "matrices must be of equal dimensions"
#        else:
#                print "Computing covariance overlap"
#
       # compute the eigenvalues and eigenvectors of matrix A
#        eigvalA, eigvecA = LA.eig(matA)          
#
       # compute the eigenvalues and eigenvectors of matrix B
#        eigvalB, eigvecB = LA.eig(matB)
#
       # initialize empty arrays to accept terms of covariance overlap measure
#        valsum = []
#        geomean = []
#        geomeansum = []
#        vecmult = []
#        vecsum = []
#	scaledcos=[]
#
 	#compute the sum of eigenvalues
#	valsum = sum(eigvalA[0:nmodes] + eigvalB[0:nmodes])
#	for i in range(nmodes):
#
      	# compute the sum of eigenvalues
      	# compute the geometric mean about two sets of eigenvalues
#      		geomean.append(math.sqrt(eigvalA[i] * eigvalA[i]))
      	# determine the inner product of two eigenvector matrices
#      		vecmult.append(eigvecA[:,i] * eigvecB[:,i])
#      		vecsum.append(sum(vecmult[i])**2)
#
        # compute the vector scaled by the geometric mean of the 
        # eigenvalues	
#	for x in range(nmodes):
#		scaledcos.append(geomean[x] * vecsum[x])
#	
#	scaledcossum = sum(scaledcos)
#
#        omegain = valsum - round((2 * scaledcossum), 20)
#        omega = 1 - (omegain/valsum)**0.5
#      
#
#	return omega

#____________________________________________________________________________
# Plotting functions 
#____________________________________________________________________________

def calcCC():
	veccosine = []
	linveccosine = []	

	for a in range(len(data)-1):
		b = a + 1
		print 'Calculating cosine content between blocks, using %s eigenmode(s):' %nmodes
		print a, b
		linveccosine.append(cosinecontent(data[a], data[b]))

	np.savetxt('time_CosCont.dat', linveccosine)	

	# plot linear cosine content
	plt.plot(range(len(linveccosine)), linveccosine, color='black')	
    	plt.grid()
    	axes = plt.gca()
	plt.ylim(0,1)
    	plt.ylabel(r'Cosine content, $\Psi_(A,B)$')
    	plt.xlabel(r'Simulation block')
    	plt.savefig('time_CosCont.pdf', bbox_inches='tight')
	
def calcCCmat():
	veccosine = []	
	for i, j in it.product(data, repeat=2):
		veccosine.append(cosinecontent(i, j))
	#np.savetxt('cosine_matrix.dat', veccosine)

	# reshape np array into matrix format
	cosinemat = np.reshape(veccosine, (len(data), len(data)))

	# convert matrix from long format to floats
	cosinemat = np.array(cosinemat, dtype=float)

	np.savetxt('cosine_matrix.dat', cosinemat)	

	plt.imshow(cosinemat, cmap='jet')
	plt.colorbar()
	#plt.clim(0,1)
	plt.savefig('cosineconent_mat.pdf')

#plt.figure()	
#calcCC()

# compute the covariance overlap for the cartesian product of all MI matrices
def calcCO():
	overlap = []
	
	last = len(data)-1	
	for a in range(len(data)-1):
		b = a + 1
		print 'Calculating cosine content between blocks, using %s eigenmode(s):' %nmodes
		print a, b
		overlap.append(covaroverlap(data[last], data[b]))
	np.savetxt('time_CovOverlap.dat', overlap)	

	# plot linear covariance overlap
	plt.plot(range(len(overlap)), overlap, color='red')	
    	plt.grid()
    	axes = plt.gca()
	plt.ylim(0,1)
    	plt.ylabel(r'Covariance overlap, $\Omega_(A,B)$')
    	plt.xlabel(r'Simulation block')
    	plt.savefig('time_CovOver.pdf', bbox_inches='tight')

def calcCOsum():
	overlap = []
	
	for a in range(len(data)-1):
		b = a + 1
		print 'Calculating cosine content between blocks, using %s eigenmode(s):' %nmodes
		print a, b
		noverlap = 1-covaroverlap(data[a], data[b])
		overlap.append(noverlap)

	np.savetxt('time_CovOverlap.dat', overlap)	

	# plot linear covariance overlap
	plt.plot(range(len(overlap)), np.cumsum(overlap)**0.5, color='red')	
    	plt.grid()
    	axes = plt.gca()
	#plt.ylim(0,1)
    	plt.ylabel(r'Cumulative covariance overlap, $\Omega_(A,B)$')
    	plt.xlabel(r'Simulation block')
    	plt.savefig('time_CovOver_cum.pdf', bbox_inches='tight')

def calcCOmat():
	overlap = []
	for x, y in it.product(data, repeat=2):
		overlap.append(covaroverlap(x, y))
		
        # reshape np array into matrix format
        covarmat = np.reshape(overlap, (len(data), len(data)))

	# convert matrix from long format into floats
	covarmat = np.array(covarmat, dtype=float)

	np.savetxt('covar_overlap_matrix.dat', covarmat, delimiter='\t')
	
	plt.imshow(covarmat, cmap='jet')
	plt.colorbar()
	#plt.clim(0,1)
	plt.savefig('covoverlap_mat.pdf')



## Calculate the first eigenvalue over time
def eigvalt():
	eigenproj = []
	
	for a in range(len(data)-1):
		print 'Calculating the first eigenvalue [%s]' %a
		eigenval, eigenvec = LA.eig(data[a])
		proj = np.dot(eigenvec[:,0], pdbdata[a])
		eigenproj.append(proj)

	# plot linear covariance overlap
	plt.plot(range(len(eigenproj)), eigenproj, color='black')	
    	plt.grid()
    	axes = plt.gca()
	#plt.ylim(0,1)
    	#plt.ylabel(r'Covariance overlap, $\Omega_(A,B)$')
    	#plt.xlabel(r'Simulation block')
    	plt.savefig('ev1_time.pdf', bbox_inches='tight')
#plt.figure()	
#eigvalt()


#____________________________________________________________________________
# Plot computed results
#____________________________________________________________________________

## SWITCH: if user selects, calculate covariance overlap
if args.CO == 'yes':
	plt.figure()
	calcCO()

## SWITCH: if user selects, calculate the cummulative sum of the covariance overlap
if args.COsum == 'yes':
	plt.figure()
	calcCOsum()



## SWITCH: if user selects matrix, compute cosine content matrix
if args.CCmat == 'yes':
	plt.figure()
	calcCCmat()


## SWITCH: if user selects matrix, compute cosine content matrix
if args.COmat == 'yes':
	plt.figure()
	calcCOmat()
