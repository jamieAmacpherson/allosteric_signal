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

gmx trjconv -f $1\
             -s $2 -n $6\
             -b $k -e $k\
	     -o $k.pdb <<EOF
3
3
EOF

gmx covar -f $k.xtc -s $k.pdb\
	  -ascii ${k}_covar.dat -o ${k}_eigval.xvg <<EOF
3
3
EOF

gmx anaeig -f $k.pdb -s $k.pdb\
	   -proj ${k}_proj.xvg -first 1 -last 1 <<EOF
3
3
EOF

rm $k.xtc eigenvec.trr average.pdb covar.log
#python $mdconvert $k.xtc -o $k.dcd


#python $covarmat -t $k.dcd -s $k.pdb -b $k 


let k=k+$3
done

#rm *.dcd

#python $overlap $5 --matrix yes
exit
