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

OUTDIR=block$k

mkdir $OUTDIR/

let b=k-$3

echo $b
echo $k
$SRCDIR/g_sa_encode -f $1\
                      -s $2\
		      -b $b -e $k\
                      -strlf $OUTDIR/$k.lf_str.out\
                      -rmsdlf $OUTDIR/$k.lf_rmsd.xvg\
                      -xpmlf $OUTDIR/$k.lf.xpm\
                      -fasta -xpm -log $OUTDIR/$k.log  

mpirun -np 6\
                $SRCDIR/g_sa_analyze -sa $OUTDIR/$k.lf_str.out\
                       -MImat $OUTDIR/$k.lf_MImat.out\
                       -eeMImat $OUTDIR/$k.lf_eeMImat.out\
                       -jHmat $OUTDIR/$k.lf_jHmat.out\
                       -nMImat $k.lf_nMImat.out\
#                       -nSample 50\
#                       -ZMImat $OUTDIR/$k.lf_ZImat.out\
#                       -meanMImat $OUTDIR/$k.lf_meanMImat.out\
#                       -stdMImat $OUTDIR/$k.lf_stdMImat.out\
#		       -pvalueMImat $OUTDIR/$k.lf_pvalueMImat.out\
                       -MImatrix -verbose 
sleep 3

let k=k+$3
done

python /home/macphej/jm.software/development/allosteric_signal/replica_ana/src/python/MI_space.py 10

Rscript /home/macphej/jm.software/development/allosteric_signal/replica_ana/src/r/MI_ana.R

exit
