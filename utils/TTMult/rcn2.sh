#!/bin/sh

# This script can be used to run the rcn2 program from the
# TTMult suite. 

if [ ! -x "$TTMULT_HOME/rcn2" ]; then
    echo "rcn2 command not found."
    exit 1
fi

if [ "$#" -eq 0 ]; then
    NAME='input'
else
    NAME="$1"
fi

# The input file for rcn2 must be added in the current folder.
if [ -f "$NAME.rcn2" ]; then
    ln -sf $NAME.rcn2 fort.10
    rcn2
    mv fort.9 $NAME.rcn2_out
    mv fort.11 $NAME.rcg.orig
    rm fort.10
    echo "rcn2 has finished successfully."
else
    echo "Could not find $NAME.rcn2 in the current folder."
    exit 1
fi

exit 0
