
cd ../03_Gridmaps/script
mkdir -p ../log
for i in `ls`; 
do
	sh $i > ../log/${i%.*}.log &
done
wait
