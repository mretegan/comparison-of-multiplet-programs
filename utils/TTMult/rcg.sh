#!/bin/sh

# This script can be used to run the ttrcg program from the
# TTMult suite. 

if [ ! -x "$TTMULT_HOME/ttrcg" ]; then
    echo "ttrcg command not found."
    exit 1
fi

if [ "$#" -eq 0 ]; then
    NAME='input'
else
    NAME="$1"
fi

# The input file $NAME.rcg has to be created using  
# $NAME.rcg.orig file as a starting point.
# See this page for more details: http://www.anorg.chem.uu.nl/CTM4XAS/tutorial_rcg.html.
if [ -f "$NAME.rcg" ]; then
    cp $TTMULT_HOME/rcg_cfp72 fort.72
    cp $TTMULT_HOME/rcg_cfp73 fort.73
    cp $TTMULT_HOME/rcg_cfp74 fort.74
    ln -sf $NAME.rcg fort.10
    ttrcg
    mv fort.9 $NAME.rcg_out
    mv fort.14 $NAME.rcg_rme
    rm fort.10 fort.72 fort.73 fort.74 
    if [ -f "FTN02" ]; then
        rm FTN02
    fi
    echo "ttrcg has finished successfully."
else
    echo "Could not find $NAME.rcg in the current folder."
    exit 1
fi

exit 0
