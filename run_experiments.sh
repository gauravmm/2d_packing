#!/bin/bash

# for FILEN in 01 02 04 10 09 08 07 06 05 03
# do
#     for PROB in {1..10}
#     do
#         timeout -k 60s 2060s julia simple_packing.jl "data/unibo/Class_${FILEN}.2bp" "$PROB" | tee "preproc-class$FILEN-prob$PROB.txt"
#     done
# done

for FILEN in 01
do
    for PROB in {1..50}
    do
        timeout -k 60s 1060s julia simple_packing.jl unibo "data/unibo/Class_${FILEN}.2bp" "$PROB" | tee "../output_nopreproc/run-class$FILEN-prob$PROB.txt"

    done
done
