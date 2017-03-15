if [ $# -lt 2 ]; then
    echo "Usage: ./run-data-with-latest-corrections.sh [5/8] [DATA_TREE] [label]"
    exit 1
fi

case $1 in
    5)
        echo "producing final results for 5.02 TeV"
        _RUN=285090
        ;;
    8)
        echo "producing final results for 8.16 TeV"
        _RUN=285832
        ;;
    *)
        echo "check sqrt(snn) energy!"
        exit 1
        ;;
esac

sed -i 's/.*define _ENERGY.*/#define _ENERGY '${1}'/' plotFinalResult.C
g++ plotFinalResult.C `root-config --cflags --libs` -O2 -Wall -o plotFinalResult

echo "producing final results..."
./plotFinalResult 12 $2 $_RUN-${3} 1 ${1}tev-latest 1 1 &> logs/$_RUN-${3}-12.log &
./plotFinalResult 13 $2 $_RUN-${3} 1 ${1}tev-latest 1 1 &> logs/$_RUN-${3}-13.log &
./plotFinalResult 23 $2 $_RUN-${3} 1 ${1}tev-latest 1 1 &> logs/$_RUN-${3}-23.log &
wait

mv correction-*-$_RUN-${3}.root correction/
./merge_tracklets $_RUN-${3} ${1}tev-latest
