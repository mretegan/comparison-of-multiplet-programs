#!/bin/sh

# This script can be used to run the ttplo program from the 
# TTMult suite. Note that the input file $NAME.plo has 
# to be created.

if [ ! -x "$TTMULT_HOME/ttplo" ]; then
    echo "ttplo command not found."
    return 1
fi

if [ "$#" -eq 0 ]; then
    NAME='input'
else
    NAME="$1"
fi

if [ -f "$NAME.plo" ]; then
    ttplo < $NAME.plo
else
    echo "Could not find $NAME.plo in the current folder."
    exit 1
fi

exit 0
