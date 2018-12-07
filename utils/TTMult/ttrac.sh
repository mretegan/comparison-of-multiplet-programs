#!/bin/sh

# This script can be used to run the ttrac program from the
# TTMult suite. RAC stands for RACER, the predecessor of RACAH.

if [ ! -x "$TTMULT_HOME/ttrac" ]; then
    echo "ttrac command not found."
    exit 1
fi

if [ "$#" -eq 0 ]; then
    NAME='input'
else
    NAME="$1"
fi

if [ -f "$NAME.rac" ]; then
    ttrac < $NAME.rac $NAME.rcg_rme $NAME.rac_out
    echo "ttrac has finished successfully."
else
    echo "Could not find $NAME.rac in the current folder."
    exit 1
fi

if [ -f rme_out.dat ]; then
    mv rme_out.dat $NAME.rac_rme
    cp $NAME.rac_rme $NAME.rac_rme.orig
    rac_converter $NAME.rac_rme.orig # real
    echo "rme_out.dat converted successfully."
else
    echo "rme_out.dat conversion failed."
    exit 1
fi

exit 0
