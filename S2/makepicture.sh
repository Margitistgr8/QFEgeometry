#!/bin/bash


make 
./affineflow
wolframscript -file makeVideo.wls 

NAME="L64"
mkdir -p ${NAME}
mv *.dat ${NAME}
mv *.png ${NAME}
