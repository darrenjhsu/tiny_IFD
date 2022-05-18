cd ..
for i in `ls`; 
do
	if [ "$i" != "script" ]; then
		cd $i
		for j in `ls *_*_*.pdb`;
		do
			echo "In $i for working on ligand $j"
			cd ${j%.*}
			rm ${j}_revert.pdb
			obabel -ipdbqt ${j}qt -opdb -O ${j%.*}_revert.pdb -d
			cd ../
		done
		cd ../
	fi
done

