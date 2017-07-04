#!/bin/bash

if [ $# -ne 6 ]
then
        echo "Incorrect number of arguments..."
        echo "Usage: gsa_sliding_window.sh <trajectory.xtc> <topology> <length of window (ps)> <final timestep (ps)> <neigen modes> <index.ndx>"
        exit 1
fi

mdconvert='/home/macphej/jm.software/development/allosteric_signal/replica_ana/src/python/mdconvert.py'
covarmat='/home/macphej/jm.software/development/allosteric_signal/replica_ana/src/python/covar_mat.py'
overlap='/home/macphej/jm.software/development/allosteric_signal/replica_ana/src/python/MI_space.py'


## Split trajectory into sliding windows 
let k=$3
let nst=$4

while [ $k -le $nst ]
do
let b=k-$3

gmx trjconv -f $1\
             -s $2 -n $6\
             -b 0 -e $k\
	     -o $k.xtc <<EOF
3
3
EOF

gmx trjconv -f $1\
             -s $2 -n $6\
             -b $k -e $k\
	     -o $k.pdb <<EOF
3
3
EOF

gmx covar -f $k.xtc -s $k.pdb\
	  -ascii ${k}_covar.dat -o ${k}_eigval.xvg -v ${k}_eigvec.trr <<EOF
3
3
EOF


gmx anaeig -f $k.pdb -s $k.pdb\
	   -v ${k}_eigvec.trr -v2 ${b}_eigvec.trr\
	   -eig ${k}_eigval.xvg -eig2 ${b}_eigval.xvg -over ${k}_overlap.xvg >& overlap_${k}.dat

#rm $k.xtc  average.pdb covar.log

sed '61q;d' overlap_${k}.dat >> covar_overlap.dat

let k=k+$3
done

rm *.trr

#python $overlap $5 --matrix yes
exit
