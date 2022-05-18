
cd ../05_Refinement/script
mkdir -p ../log
for i in `ls prepareComplex{?,??}.sh`; 
do
	sh $i > ../log/${i%.*}.log &
done

wait

for i in `ls prepareComplex_openmm*.sh`; 
do
	sh $i > ../log/${i%.*}.log &
done

wait
# Also run python PlanEM.py

#sh planEM.sh

