./ana2.sh 12 $1  "Data-$2" 1 $3 10000000000 $4 $5
mv correction.root plot/root/result-12-$2-$3.root
./ana2.sh 23 $1  "Data-$2" 1 $3 10000000000 $4 $5
mv correction.root plot/root/result-23-$2-$3.root
./ana2.sh 13 $1  "Data-$2" 1 $3 10000000000 $4 $5
mv correction.root plot/root/result-13-$2-$3.root

./merge.sh Data-$2
root -l -q -b "makeMergedPlot.C+(\"$2-$3\",0,3.8,($4+$5)/2)"
