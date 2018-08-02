#!/bin/bash

if [ $# -ne 3 ]
then
        echo "Incorrect number of arguments..."
        echo "Usage: gsa_encode.sh <trajectory.xtc> <topology.tpr> <output prefix>"
        exit 1
fi

set -e

source /usr/local/gromacs/bin/GMXRC.bash

SRCDIR=/home/macphej/jm.software/apps/gsatools-4.5.x-1.00/src

OUTPREFIX=$3
OUTDIR=.
###################
###################
$SRCDIR/g_sa_encode -f $2\
                      -s $2\
                      -strlf $OUTDIR/$OUTPREFIX.lf_str.out\
                      -rmsdlf $OUTDIR/$OUTPREFIX.lf_rmsd.xvg\
                      -xpmlf $OUTDIR/$OUTPREFIX.lf.xpm\
                      -fasta -xpm -log $OUTDIR/$OUTPREFIX.log & 
