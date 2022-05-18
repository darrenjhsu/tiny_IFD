
cd ../04_Docking/script
mkdir -p ../log
for i in `ls dock*.sh`; 
do
  echo $i
  sleep 3
  sh $i > ../log/${i%.*}.log
done
