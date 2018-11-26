#!/bin/sh

# This script can be used to run the RCN2 program from the
# TTMult suite. 

if [ ! -x "$TTMULT_HOME/rcn2" ]; then
    echo "rcn2 command not found."
    exit 1
fi

if [ "$#" -eq 0 ]; then
    echo 'You need to provide the base name as the first argument.'
    exit 1
fi

NAME=$1

# The input file for RCN2 must be added in the current folder.
if [ -f "$NAME.rcn2" ]; then
    ln -sf $NAME.rcn2 fort.10
    rcn2
    mv fort.9 $NAME.rcn2_out
    mv fort.11 $NAME.rcg.orig
    rm fort.10
    echo "RCN2 has finished successfully."
    exit 0
else
    echo "Could not find $NAME.rcn2 in the current folder."
    exit 1
fi
