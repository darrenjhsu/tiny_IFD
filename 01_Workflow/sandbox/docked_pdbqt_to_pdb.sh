
cd ..
for i in `ls`; 
do
	if [ "$i" != "dock.sh" ]; then
		cd $i
		for j in `ls *_*_*.pdb`;
		do
			echo "In $i for working on ligand $j"
			cd ${j%.*}
			mkdir docked_pdb
			for k in `ls docked/`;
			do
				obabel -ipdbqt docked/$k -opdb -O docked_pdb/${k%.*}.pdb -d
			done
			cd ../
		done
		cd ../
	fi
done

