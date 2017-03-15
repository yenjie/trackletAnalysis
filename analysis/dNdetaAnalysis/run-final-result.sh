if [ $# -lt 2 ]; then
    echo "Usage: ./run-final-result.sh [5/8] [label]"
    exit 1
fi

case $1 in
    5)
        echo "producing final results for 5.02 TeV"

        _RUN=285090
        _MC_TTREE="/data/biran/CMSSW_8_0_23/src/trackletAnalysis/analysis/trackletProduction/TrackletTree-EPOS-5TeV-LATEST.root"
        _DATA_TTREE="/data/biran/CMSSW_8_0_23/src/trackletAnalysis/analysis/trackletProduction/TrackletTree-285090-ALL.root"
        ;;
    8)
        echo "producing final results for 8.16 TeV"

        _RUN=285832
        _MC_TTREE="/data/biran/CMSSW_8_0_23/src/trackletAnalysis/analysis/trackletProduction/TrackletTree-EPOS-8TeV-LATEST.root"
        _DATA_TTREE="/data/biran/CMSSW_8_0_23/src/trackletAnalysis/analysis/trackletProduction/TrackletTree-285832-3M.root"
        ;;
    *)
        echo "check sqrt(snn) energy!"
        exit 1
        ;;
esac

sed -i 's/.*define _ENERGY.*/#define _ENERGY '${1}'/' plotFinalResult.C
g++ plotFinalResult.C `root-config --cflags --libs` -O2 -Wall -o plotFinalResult || exit 1

cp ./rootfiles/accep/${1}tev/100/acceptance-*.root ./correction/

echo "deriving correction factors..."
./plotFinalResult 12 $_MC_TTREE ${1}tev-${2} 0 ${1}tev-${2} 1 1 &> logs/${1}tev-${2}-12.log &
./plotFinalResult 13 $_MC_TTREE ${1}tev-${2} 0 ${1}tev-${2} 1 1 &> logs/${1}tev-${2}-13.log &
./plotFinalResult 23 $_MC_TTREE ${1}tev-${2} 0 ${1}tev-${2} 1 1 &> logs/${1}tev-${2}-23.log &
wait

echo "producing final results..."
mv correction-*-${1}tev-${2}.root correction/
./plotFinalResult 12 $_DATA_TTREE $_RUN-${2} 1 ${1}tev-${2} 1 1 &> logs/$_RUN-${2}-12.log &
./plotFinalResult 13 $_DATA_TTREE $_RUN-${2} 1 ${1}tev-${2} 1 1 &> logs/$_RUN-${2}-13.log &
./plotFinalResult 23 $_DATA_TTREE $_RUN-${2} 1 ${1}tev-${2} 1 1 &> logs/$_RUN-${2}-23.log &
wait

mv correction-*-$_RUN-${2}.root correction/
./merge_tracklets $_RUN-${2} ${1}tev-${2}
