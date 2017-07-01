#____________________________________________________________________________
# Calculate the positional covariance matrix for an MD simulation
#
# Copyright 2016 King's College London and the Authors
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
from prody import *
from pylab import *
import MDAnalysis as mda
from MDAnalysis.analysis.align import *
import argparse
import os.path
import sys
import numpy as np
import subprocess
import itertools as it
from numpy import linalg as LA
import math as math
import array as ar
import matplotlib.pyplot as plt


#____________________________________________________________________________
# Parse commandline arguments
#____________________________________________________________________________
# check if input file exists
# there are two inputs: trajectory (.dcd) and topology (.pdb)
# if either of those inputs are not supplied, or if the user doesn't invoke the
# help flag the program will display an error message.
def is_valid_file(arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist! Use the --help flag for input options." % arg)
    else:
        return arg

# command line argument parser
parser = argparse.ArgumentParser(description='Calculate entropy of MD trajectory')

# the first argument is the trajectory file (.dcd) supplied after the -t flag
# the trajectory file is saved as an object with the variable args.dcdfile
parser.add_argument("-t", dest="dcdfile", required=True,
                    help="Free trajectory file (format: .dcd)",
                    type=lambda x: is_valid_file(x))

# the second argument is the topology file (.pdb) supplied after the -s flag
# this is saved an an obect with the variable args.pdbfile
parser.add_argument("-s", dest="pdbfile", required=True,
                    help="Free structure file (format: .pdb)",
                    type=lambda x: is_valid_file(x))

parser.add_argument("-b", dest="block", required=True)

# the arguments are parsed 
args = parser.parse_args()


#____________________________________________________________________________
# remove rotational-translational motions from trajectory
#____________________________________________________________________________
# to calculate the configurational entropy, we first remove the rotational-translational
# motions from the trajectory. This is done by fitting the trajectory to the reference
# structure (ie. the pdbfile). 

def rmrt(topology, trajectory):
    # define the reference structure as the topology file
    ref = mda.Universe(topology)
#    
    # define the trajectory as the .dcd file and link it to the topology symbolically
    traj = mda.Universe(topology, trajectory)
#    
    # fit the trajectory to the topology, removing rotational-translational motions
    # and save the resulting trajectory.
    #rms_fit_trj(traj, ref, filename='rmsfit_traj.dcd')

rmrt(args.pdbfile, args.dcdfile)


#____________________________________________________________________________
# covariance matrix
#____________________________________________________________________________
def covar(topology, trajectory, block):
    struct = parsePDB(topology)
    traj = Trajectory(trajectory)
    traj.link(struct)
    traj.setCoords(struct)
    traj.setAtoms(struct.calpha)
    ensemble = EDA('trajectory')
    ensemble.buildCovariance( traj )
    mat = ensemble.getCovariance()
    np.savetxt('%s.out' % block, mat)

covar(args.pdbfile, args.dcdfile, args.block)

