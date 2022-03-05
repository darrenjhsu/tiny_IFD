
cd ../05_Refinement/script
mkdir -p ../log
for i in `ls antechamber*`; 
do
	sh $i > ../log/${i%.*}.log &
done

wait

# Also run python PlanEM.py

#sh planEM.sh

