#!/bin/sh

# This script can be used to run the ttban program from the 
# TTMult suite. Note that the input file $NAME.ban has 
# to be created.

if [ ! -x "$TTMULT_HOME/ttban" ]; then
    echo "ttban command not found."
    return 1
fi

if [ "$#" -eq 0 ]; then
    NAME='input'
else
    NAME="$1"
fi


if [ -f $NAME.ban ]; then
    cp $NAME.ban fort.50
    ln -sf $NAME.rcg_rme FTN14
    ln -sf $NAME.rac_rme FTN15
    ttban
    mv fort.44 $NAME.ban_out
    rm FTN14 FTN15 FTN98 FTN99 
    rm -fr fort.*
else
    echo "Could not find $NAME.ban in the current folder."
    exit 1
fi

exit 0
