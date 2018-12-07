#!/bin/sh

# This script can be used to run the rcn program from the 
# TTMult suite. Note that the input file $NAME.rcn has 
# to be created.

if [ ! -x "$TTMULT_HOME/rcn" ]; then
    echo "rcn command not found."
    return 1
fi

if [ "$#" -eq 0 ]; then
    NAME='input'
else
    NAME="$1"
fi

if [ -f "$NAME.rcn" ]; then
    ln -sf $NAME.rcn fort.10
    # Note that rcn is an alias to rcn31.
    rcn 
    mv fort.9 $NAME.rcn_out
    rm fort.10
    echo "rcn has finished successfully."
else
    echo "Could not find $NAME.rcn in the current folder."
    exit 1
fi

exit 0
