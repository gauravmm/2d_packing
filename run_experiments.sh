#!/bin/bash

for FILEN in 01 02 04 10 09 08 07 06 05 03
do
    for PROB in {1..10}
    do
        echo julia simple_packing.jl "data/unibo/Class_${FILEN}.2bp" "$PROB" | tee "preproc-class$FILEN-prob$PROB.txt"
    done
done
