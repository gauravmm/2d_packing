#!/bin/bash

for FILEN in 09
do
    for PROB in {11..50}
    do
        timeout -k 60s 2060s julia simple_packing.jl unibo "data/unibo/Class_${FILEN}.2bp" "$PROB" | tee "../output_preproc/run-class$FILEN-prob$PROB.txt"
    done
done
