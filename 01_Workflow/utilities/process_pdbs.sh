#!/bin/bash

for i in `ls *.pdb`;
do
  sh ../antechamber_pipeline-master/pipeline.sh $i
done
