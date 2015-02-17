#echo all energies

# first and second layers
./ana.sh 12 $1  "900GeV MC 12" 0 
mv correction.root correction/correction-12-$2.root

# second and third layers
./ana.sh 23 $1  "900GeV MC 23" 0 
mv correction.root correction/correction-23-$2.root

# first and third layers
./ana.sh 13 $1  "900GeV MC 13" 0
mv correction.root correction/correction-13-$2.root


# Get Trigger Correction
#./trig.sh $1 $2
