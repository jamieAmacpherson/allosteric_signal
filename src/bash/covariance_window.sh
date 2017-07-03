#!/bin/bash

let b=k-$3

echo $b
echo $k

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



echo $b
echo $k
gmx trjconv -f $1\
             -s $2 -n $6\
             -b $b -e $k\
	     -o $b.xtc <<EOF
3
3
EOF


python $mdconvert $b.xtc -o $b.dcd

rm $b.xtc 

python $covarmat -t $b.dcd -s topol.pdb -b $b 


let k=k+$3
done

rm *.dcd

python $overlap $5 --matrix yes
exit

