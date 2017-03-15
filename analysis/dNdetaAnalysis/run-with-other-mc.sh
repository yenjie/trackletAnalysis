if [ $# -lt 3 ]; then
    echo "Usage: ./run-with-other-mc.sh [5/8] [label] [MC TTree]"
    exit 1
fi

if [ $1 -eq 5 ]; then
    echo "producing final results for 5.02 TeV"

    _RUN=285090
    _DATA_TTREE="/data/biran/CMSSW_8_0_23/src/trackletAnalysis/analysis/trackletProduction/TrackletTree-285090-ALL.root"
elif [ $1 -eq 8 ]; then
    echo "producing final results for 8.16 TeV"

    _RUN=285832
    _DATA_TTREE="/data/biran/CMSSW_8_0_23/src/trackletAnalysis/analysis/trackletProduction/TrackletTree-285832-3M.root"
else
    echo "check sqrt(snn) energy!"
fi

sed -i 's/.*define _ENERGY.*/#define _ENERGY '${1}'/' plotFinalResult.C
g++ plotFinalResult.C `root-config --cflags --libs` -O2 -Wall -o plotFinalResult

echo "deriving correction factors..."
./plotFinalResult 12 ${3} ${1}tev-${2} 0 ${1}tev-${2} 1 1 &> logs/${1}tev-${2}-12.log &
./plotFinalResult 13 ${3} ${1}tev-${2} 0 ${1}tev-${2} 1 1 &> logs/${1}tev-${2}-13.log &
./plotFinalResult 23 ${3} ${1}tev-${2} 0 ${1}tev-${2} 1 1 &> logs/${1}tev-${2}-23.log &
wait

echo "producing final results..."
mv correction-*-${1}tev-${2}.root correction/
./plotFinalResult 12 $_DATA_TTREE $_RUN-${2} 1 ${1}tev-${2} 1 1 &> logs/$_RUN-${2}-12.log &
./plotFinalResult 13 $_DATA_TTREE $_RUN-${2} 1 ${1}tev-${2} 1 1 &> logs/$_RUN-${2}-13.log &
./plotFinalResult 23 $_DATA_TTREE $_RUN-${2} 1 ${1}tev-${2} 1 1 &> logs/$_RUN-${2}-23.log &
wait

mv correction-*-$_RUN-${2}.root correction/
./merge_tracklets $_RUN-${2} ${1}tev-${2}
