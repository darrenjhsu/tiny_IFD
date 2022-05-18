#!/bin/bash

cd ../
for i in `ls -d *-*-*`
do
  echo $i `cat ${i}/Structure/prmtop/${i}.prmtop | head -7 | tail -1 | awk '{printf $1}'`
done
