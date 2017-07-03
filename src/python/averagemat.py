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

parser = argparse.ArgumentParser(description='Calculate the covariance overlap and cosine content of Mutual information matrices.')


# the first argument is the trajectory file (.dcd) supplied after the -t flag
# the trajectory file is saved as an object with the variable args.dcdfile
parser.add_argument("-a1", dest="mata1", required=True,
                    help="First Hermitian matrix",
                    type=lambda x: is_valid_file(x))

parser.add_argument("-a2", dest="mata2", required=True,
                    help="Second Hermitian matrix",
                    type=lambda x: is_valid_file(x))

parser.add_argument("-a3", dest="mata3", required=True,
                    help="Second Hermitian matrix",
                    type=lambda x: is_valid_file(x))

parser.add_argument("-b1", dest="matb1", required=True,
                    help="Second Hermitian matrix",
                    type=lambda x: is_valid_file(x))

parser.add_argument("-b2", dest="matb2", required=True,
                    help="Second Hermitian matrix",
                    type=lambda x: is_valid_file(x))

parser.add_argument("-b3", dest="matb3", required=True,
                    help="Second Hermitian matrix",
                    type=lambda x: is_valid_file(x))

args = parser.parse_args()


#____________________________________________________________________________
# Compute the average MI matrices and take the difference matrix
#____________________________________________________________________________

def avdiff(mat1, mat2, mat3, mat4, mat5, mat6):
	mat1 = np.loadtxt(mat1)
	mat2 = np.loadtxt(mat2)
	mat3 = np.loadtxt(mat3)
	
	mat4 = np.loadtxt(mat4)
	mat5 = np.loadtxt(mat5)
	mat6 = np.loadtxt(mat6)

	av1 = (mat1 + mat2 + mat3)/3
	av2 = (mat4 + mat5 + mat6)/3

	np.savetxt('mat1av.mat', av1)
	np.savetxt('mat2av.mat', av2)
#	if np.shape(av1) =! np.shape(av2):
#		av1row, av1col = np.shape(av1)
#		av2row, av2col = np.shape(av2)
#
#		difrow = av1row - av2row
#		difcol = av1col - av2col
#			
#		if difrow < 0:
#			difrow = difrow * -1
#
#		if difcol < 0:
#			difcol = difcol * -1
#
#		if av1row > av2row:
#			av1 = np.delete(av1, av1row, 0)
#		elif av2row > av1row:
#			av2 = np.
		
	diffmat = av1 - av2
	np.savetxt('difference.mat', diffmat)


avdiff(args.mata1, args.mata2, args.mata3, args.matb1, args.matb2, args.matb3)



