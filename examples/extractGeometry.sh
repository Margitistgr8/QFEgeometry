#!/bin/bash

exec=/Users/jinyunlin/QFE-Research/QFEgeometry/bin/S2ReadGeometry
q=5
LList=()
#for i in 4 8 10 16 20 24 32 36 48 52 64 72 84 96
for i in 8 10 24 36 52 72 84
do 
Data_Dir=/Users/jinyunlin/QFE-Research/QFEgeometry/data/q${q}/xivec_${i}.dat
echo ${Data_Dir}
${exec} ${i} ${q} ${Data_Dir} 
done 