#!/bin/bash

IFILE='instruction_xfoil.txt'

NACA_LIST=("0006" "0012")

for NACA_CASE in ${NACA_LIST[@]}
do 
    for ALPHA in {1..2}
    do
        echo "NACA $NACA_CASE" > $IFILE
        echo "OPER" >> $IFILE
        echo "ALFA $ALPHA" >> $IFILE
        echo "CPWR naca_"$NACA_CASE"_"$ALPHA"_cp.txt" >> $IFILE 
        echo '  ' >> $IFILE
        echo "QUIT" >> $IFILE
        xfoil < $IFILE
    done
    echo ''
done 



    