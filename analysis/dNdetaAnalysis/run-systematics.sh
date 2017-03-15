if [ $# -lt 1 ]; then
    echo "Usage: ./run-systematics.sh [5/8]"
    exit 1
fi

case $1 in
    5)
        echo "producing systematics for 5.02 TeV"

        _RUN=285090
        _MC_TTREE="/data/biran/CMSSW_8_0_23/src/trackletAnalysis/analysis/trackletProduction/TrackletTree-EPOS-5TeV-LATEST.root"
        _MC_NOISE="/data/biran/CMSSW_8_0_23/src/trackletAnalysis/analysis/trackletProduction/TrackletTree-EPOS-5TeV-NOISE.root"
        _DATA_TTREE="/data/biran/CMSSW_8_0_23/src/trackletAnalysis/analysis/trackletProduction/TrackletTree-285090-ALL.root"
        _SPLIT_TTREE="/data/biran/CMSSW_8_0_23/src/trackletAnalysis/analysis/trackletProduction/TrackletTree-285090-SPLIT.root"
        _DROP_TTREE="/data/biran/CMSSW_8_0_23/src/trackletAnalysis/analysis/trackletProduction/TrackletTree-285090-DROP.root"
        _SMEAR_TTREE="/data/biran/CMSSW_8_0_23/src/trackletAnalysis/analysis/trackletProduction/TrackletTree-285090-SMEAR.root"
        ;;
    8)
        echo "producing final results for 8.16 TeV"

        _RUN=285832
        _MC_TTREE="/data/biran/CMSSW_8_0_23/src/trackletAnalysis/analysis/trackletProduction/TrackletTree-EPOS-8TeV-LATEST.root"
        _MC_NOISE="/data/biran/CMSSW_8_0_23/src/trackletAnalysis/analysis/trackletProduction/TrackletTree-EPOS-8TeV-NOISE.root"
        _DATA_TTREE="/data/biran/CMSSW_8_0_23/src/trackletAnalysis/analysis/trackletProduction/TrackletTree-285832-3M.root"
        _SPLIT_TTREE="/data/biran/CMSSW_8_0_23/src/trackletAnalysis/analysis/trackletProduction/TrackletTree-285832-SPLIT.root"
        _DROP_TTREE="/data/biran/CMSSW_8_0_23/src/trackletAnalysis/analysis/trackletProduction/TrackletTree-285832-DROP.root"
        _SMEAR_TTREE="/data/biran/CMSSW_8_0_23/src/trackletAnalysis/analysis/trackletProduction/TrackletTree-285832-SMEAR.root"
        ;;
    *)
        echo "check sqrt(snn) energy!"
        exit 1
esac

sed -i 's/.*define _ENERGY.*/#define _ENERGY '${1}'/' plotFinalResult.C
g++ plotFinalResult.C `root-config --cflags --libs` -O2 -Wall -o plotFinalResult

cp ./rootfiles/accep/${1}tev/100/acceptance-*.root ./correction/

echo "running noise systematics..."
./run-with-other-mc.sh ${1} noise-sys $_MC_NOISE
echo "running pixel splitting systematics..."
./run-data-with-latest-corrections.sh ${1} $_SPLIT_TTREE split-sys
echo "running pixel inefficiency systematics..."
./run-data-with-latest-corrections.sh ${1} $_DROP_TTREE drop-sys
echo "running misalignment systematics..."
./run-data-with-latest-corrections.sh ${1} $_SMEAR_TTREE smear-sys

echo "running mult parametrisation (nhit1_cut) systematics..."
./plotFinalResult 12 $_MC_TTREE ${1}tev-mult-1-sys 0 ${1}tev-mult-1-sys 1 1 1 1 1 &> logs/${1}tev-mult-1-sys-12.log &
./plotFinalResult 13 $_MC_TTREE ${1}tev-mult-1-sys 0 ${1}tev-mult-1-sys 1 1 1 1 1 &> logs/${1}tev-mult-1-sys-13.log &
./plotFinalResult 23 $_MC_TTREE ${1}tev-mult-1-sys 0 ${1}tev-mult-1-sys 1 1 1 1 1 &> logs/${1}tev-mult-1-sys-23.log &
wait

mv correction-*-${1}tev-mult-1-sys.root correction/
./plotFinalResult 12 $_DATA_TTREE $_RUN-mult-1-sys 1 ${1}tev-mult-1-sys 1 1 1 1 1 &> logs/$_RUN-mult-1-sys-12.log &
./plotFinalResult 13 $_DATA_TTREE $_RUN-mult-1-sys 1 ${1}tev-mult-1-sys 1 1 1 1 1 &> logs/$_RUN-mult-1-sys-13.log &
./plotFinalResult 23 $_DATA_TTREE $_RUN-mult-1-sys 1 ${1}tev-mult-1-sys 1 1 1 1 1 &> logs/$_RUN-mult-1-sys-23.log &
wait

mv correction-*-$_RUN-mult-1-sys.root correction/
./merge_tracklets $_RUN-mult-1-sys ${1}tev-mult-1-sys

echo "running mult parametrisation (nTracklets) systematics..."
./plotFinalResult 12 $_MC_TTREE ${1}tev-mult-2-sys 0 ${1}tev-mult-2-sys 1 1 1 1 2 &> logs/${1}tev-mult-2-sys-12.log &
./plotFinalResult 13 $_MC_TTREE ${1}tev-mult-2-sys 0 ${1}tev-mult-2-sys 1 1 1 1 2 &> logs/${1}tev-mult-2-sys-13.log &
./plotFinalResult 23 $_MC_TTREE ${1}tev-mult-2-sys 0 ${1}tev-mult-2-sys 1 1 1 1 2 &> logs/${1}tev-mult-2-sys-23.log &
wait

mv correction-*-${1}tev-mult-2-sys.root correction/
./plotFinalResult 12 $_DATA_TTREE $_RUN-mult-2-sys 1 ${1}tev-mult-2-sys 1 1 1 1 2 &> logs/$_RUN-mult-2-sys-12.log &
./plotFinalResult 13 $_DATA_TTREE $_RUN-mult-2-sys 1 ${1}tev-mult-2-sys 1 1 1 1 2 &> logs/$_RUN-mult-2-sys-13.log &
./plotFinalResult 23 $_DATA_TTREE $_RUN-mult-2-sys 1 ${1}tev-mult-2-sys 1 1 1 1 2 &> logs/$_RUN-mult-2-sys-23.log &
wait

mv correction-*-$_RUN-mult-2-sys.root correction/
./merge_tracklets $_RUN-mult-2-sys ${1}tev-mult-2-sys

echo "running dphi region systematics..."
./plotFinalResult 12 $_MC_TTREE ${1}tev-dphi-sys 0 ${1}tev-dphi-sys 1 1 1 1 0 1 &> logs/${1}tev-dphi-sys-12.log &
./plotFinalResult 13 $_MC_TTREE ${1}tev-dphi-sys 0 ${1}tev-dphi-sys 1 1 1 1 0 1 &> logs/${1}tev-dphi-sys-13.log &
./plotFinalResult 23 $_MC_TTREE ${1}tev-dphi-sys 0 ${1}tev-dphi-sys 1 1 1 1 0 1 &> logs/${1}tev-dphi-sys-23.log &
wait

mv correction-*-${1}tev-dphi-sys.root correction/
./plotFinalResult 12 $_DATA_TTREE $_RUN-dphi-sys 1 ${1}tev-dphi-sys 1 1 1 1 0 1 &> logs/$_RUN-dphi-sys-12.log &
./plotFinalResult 13 $_DATA_TTREE $_RUN-dphi-sys 1 ${1}tev-dphi-sys 1 1 1 1 0 1 &> logs/$_RUN-dphi-sys-13.log &
./plotFinalResult 23 $_DATA_TTREE $_RUN-dphi-sys 1 ${1}tev-dphi-sys 1 1 1 1 0 1 &> logs/$_RUN-dphi-sys-23.log &
wait

mv correction-*-$_RUN-dphi-sys.root correction/
./merge_tracklets $_RUN-dphi-sys ${1}tev-dphi-sys
