
cd ..
for i in `ls` ;
do
	if [ "$i" != "script" ]; then
		cd $i
                IFS='_' read -ra PAIR <<< "$i"
		ref=${PAIR[0]}
		target=${PAIR[1]}
		
		for j in `ls *_*_*.pdb`;
		do
			rm -r ${j%.*}
			echo -e "In receptor $i docking $j from $ref to $target ...\n\n\n"
			mkdir -p ${j%.*}
			obabel -ipdb $j -opdbqt -O ${j%.*}/${j}qt -p --partialcharge eem
			cd ${j%.*}
			../../../../../AutoDock-GPU/bin/autodock_gpu_64wi -nrun 50 -autostop 1 -nev 3000000 -ffile ../${target}.fld -lfile ${j}qt
			cd ..
		done
		cd ../
	fi
done
