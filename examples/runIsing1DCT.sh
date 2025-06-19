#!/bin/bash

EXEC=/Users/jinyunlin/QFE-Research/QFEgeometry/bin/Ising1DCT

Nx=12
invT=12.0
k1=1.0
n_iter=50000
n_therm=1000
n_skip=10
baseseed=17142
i=0


for tc in 0.005 0.025 1.8 2.2
#for tc in 0.4 0.8 1.0 1.2 1.4
do
    seed=$((baseseed+i))
    $EXEC -x $Nx \
          -T $invT -a $k1 \
          -b $tc -h $n_iter -t $n_therm -s $n_skip \
          -S $seed -d $OUTPUT_DIR
done