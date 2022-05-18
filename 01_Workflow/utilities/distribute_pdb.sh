#!/bin/bash

cd SDF_match_PDB
counter=0
overall_counter=0
for i in `ls *_0.pdb`;
do
  cp $i ../pipeline$counter/${i%_*}.pdb
  counter=$(( counter + 1 ))
  counter=$(( counter % 8 ))
  overall_counter=$(( overall_counter + 1 ))
done
echo $overall_counter
