#!/bin/sh

#How many files for marker panels will there be?
Files=4


#Create a vector of panel names
PanelNames=$(for i in $(seq 1 $Files); do echo "MarkerPanel$i.csv"; done)

for n in $(seq 0 $(($Files-1)))
do
./R
library("PseudoBabies")
setwd("~/Desktop/Chapter3-Parenetage/RpackageDev_Simulations/AukeCreekSims1")
Rn.Colony(ColonyDir="/Users/patdbarry/Desktop/GeneticsSoftware/colony2",nSim=10, Geno_Error=T,ErrorVals=3,Miss_data=F,Markerfile = $PanelNames[n],ShowProgress=F)
pids="$pids $!"
done
wait $pids
