#!/bin/bash

cd ..
for i in `ls`;
do
	if [ "$i" != "script" ]; then
		cd $i
                IFS='_' read -ra PAIR <<< "$i"
		ref=${PAIR[0]}
		target=${PAIR[1]}
		echo "In folder $i processing $target"
		# Use reduce to add hydrogen
		../../../../reduce/bin/bin/reduce -NOFLIP ${target}_apo.pdb > ${target}_h.pdb
		obabel -ipdb ${target}_h.pdb -opdbqt -O monomer_${target}.pdbqt --partialcharge gasteiger -xr -xp -xc
		../../../../AutoGrid4/autodocksuite-4.0.1/bin/autogrid4 -p ${target}.gpf
		cd ..
	fi
done
