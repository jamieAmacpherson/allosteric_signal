#!/bin/bash

if [ $# -ne 4 ]
then
        echo "Incorrect number of arguments..."
        echo "Usage: gsa_sliding_window.sh <trajectory.xtc> <topology> <length of window (ps)> <final timestep (ps)>"
        exit 1
fi


source /usr/local/gromacs/bin/GMXRC.bash

SRCDIR=/home/macphej/jm.software/apps/gsatools-4.5.x-1.00/src


let k=$3
let nst=$4

while [ $k -le $nst ]
do


let b=k-$3

echo $b
echo $k
gmx covar -f $1\
             -s $2\
             -b $b -e $k\
	     -ascii $b.covar.dat <<EOF
3
3
EOF
rm average.pdb eigenval.xvg eigenvec.trr covar.log


let k=k+$3
done
exit
