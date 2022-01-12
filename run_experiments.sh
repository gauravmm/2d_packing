#!/bin/bash

for FILEN in 02
do
    for PROB in {1..50}
    do
        timeout -k 60s 2060s julia simple_packing.jl unibo "data/unibo/Class_${FILEN}.2bp" "$PROB" | tee "../output_cplex/run-class$FILEN-prob$PROB.txt"
    done
done
