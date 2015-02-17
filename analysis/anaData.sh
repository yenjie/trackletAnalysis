./ana.sh 12 $1  "Data-$2" 1 $3 
mv correction.root plot/root/result-12-$2-$3.root
./ana.sh 23 $1  "Data-$2" 1 $3 
mv correction.root plot/root/result-23-$2-$3.root
./ana.sh 13 $1  "Data-$2" 1 $3 
mv correction.root plot/root/result-13-$2-$3.root

./merge.sh Data-$2
root -l -q "makeMergedPlot.C+(\"$2-$3\")"
