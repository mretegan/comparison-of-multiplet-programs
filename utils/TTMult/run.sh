#!/bin/bash

if [[ $# -eq 0 ]]; then
    echo 'You need to provide the base name as the first argument.'
    exit 0
fi

NAME=$1

# RCN
# Create the input file $NAME.rcn 
if [[ -f $NAME.rcn ]]; then
    ln -sf $NAME.rcn fort.10
    rcn
    mv fort.9 $NAME.rcn_out
    rm fort.10
fi

# RCN2
# Copy the appropriate input file.
if [[ -f $NAME.rcn2 ]]; then
    ln -sf $NAME.rcn2 fort.10
    rcn2
    mv fort.9 $NAME.rcn2_out
    mv fort.11 $NAME.rcg.orig
    rm fort.10
fi

# RCG
# Create the input file $NAME.rcg or edit the $NAME.rcg.orig file.
cp $TTMULT_HOME/rcg_cfp72 fort.72
cp $TTMULT_HOME/rcg_cfp73 fort.73
cp $TTMULT_HOME/rcg_cfp74 fort.74
ln -sf $NAME.rcg fort.10
ttrcg
mv fort.9 $NAME.rcg_out
mv fort.14 $NAME.rcg_rme
rm fort.10 fort.72 fort.73 fort.74 
rm FTN02
 
# Racah/Racer
# Create the input file $NAME.rac.
ttrac < $NAME.rac $NAME.rcg_rme $NAME.rac_out
if [[ -f rme_out.dat ]]; then
    mv rme_out.dat $NAME.rac_rme
    cp $NAME.rac_rme $NAME.rac_rme.orig
    rac_converter $NAME.rac_rme.orig real
fi

if [[ -f $NAME.ban ]]; then
    cp $NAME.ban fort.50
    ln -sf $NAME.rcg_rme FTN14
    ln -sf $NAME.rac_rme FTN15
    ttban
    mv fort.44 $NAME.ban_out
    rm FTN14 FTN15 FTN98 FTN99 
    rm -fr fort.*
fi

# Plo
ttplo < $NAME.plo
