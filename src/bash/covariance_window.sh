#!/bin/bash

if [ $# -ne 6 ]
then
        echo "Incorrect number of arguments..."
        echo "Usage: gsa_sliding_window.sh <trajectory.xtc> <topology> <length of window (ps)> <final timestep (ps)> <neigen modes> <index.ndx>"
        exit 1
fi

mdconvert='/home/macphej/jm.software/development/allosteric_signal/src/python/mdconvert.py'
covarmat='/home/macphej/jm.software/development/allosteric_signal/src/python/covar_mat.py'
overlap='/home/macphej/jm.software/development/allosteric_signal/src/python/MI_space.py'

gmx trjconv -f $1\
             -s $2 -n $6\
             -b 0 -e 0\
	     -o topol.pdb <<EOF
3
3
EOF

## Split trajectory into blocks
let k=$3
let nst=$4

while [ $k -le $nst ]
do

gmx trjconv -f $1\
             -s $2 -n $6\
             -b 0 -e $k\
	     -o $k.xtc <<EOF
3
3
EOF


python $mdconvert $k.xtc -o $k.dcd


python $covarmat -t $k.dcd -s topol.pdb -b $k 

rm *.dcd *.xtc

let k=k+$3
done


python $overlap $5 --CO yes
exit
