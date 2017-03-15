#!/bin/bash

if [ $# -lt 3 ]; then
    echo "Usage: ./produce_analysis_samples [5/8] [0/1/2] [label]"
    echo "       0: ALL, 1: DATA, 2: MC"
    exit 1
fi

case $1 in
    5)
        echo "producing analysis tracklettree for 5.02 TeV"

        _RUN=285090
        _MC_PTREE="/data/biran/pixeltrees/PixelTree-EPOS-5TeV-2017.root"
        _DATA_PTREE="/data/biran/pixeltrees/2016pPb/Express/PixelTree-285090.root"
        _MC_SPLIT=(0 80000 160000 240000 320000 400000 480000)
        _DATA_SPLIT=(0 320000 640000 960000 1280000 1600000 1920000)
        ;;
    8)
        echo "producing analysis tracklettree for 8.16 TeV"

        _RUN=285832
        _MC_PTREE="/data/biran/pixeltrees/PixelTree-EPOS-8TeV-2017.root"
        _DATA_PTREE="/data/biran/pixeltrees/2016pPb/Express/PixelTree-285832.root"
        _MC_SPLIT=(0 135000 270000 405000 540000 675000 810000)
        _DATA_SPLIT=(3000000 3500000 4000000 4500000 5000000 5500000 6000000)
        ;;
    *)
        echo "check sqrt(snn) energy!"
        exit 1
        ;;
esac

if [ $2 -ne 1 ]; then
    ./configure.sh $1 1

    echo "compilation test"
    root -l -b -q 'analyze_trackletTree.C++("'$_MC_PTREE'", "TrackletTree-COMPILE-DUMMY.root", 0, 1)'
    rm TrackletTree-COMPILE-DUMMY.root

    echo "producing mc tracklettree..."
    root -l -b -q 'analyze_trackletTree.C+("'$_MC_PTREE'", "TrackletTree-'$1'TeV-'$3'-0.root", '${_MC_SPLIT[0]}', '${_MC_SPLIT[1]}')' &> logs/${1}tev_mc_0.log && root -l -b -q 'analyze_trackletTree.C+("'$_MC_PTREE'", "TrackletTree-'$1'TeV-'$3'-1.root", '${_MC_SPLIT[1]}', '${_MC_SPLIT[2]}')' &> logs/${1}tev_mc_1.log &

    root -l -b -q 'analyze_trackletTree.C+("'$_MC_PTREE'", "TrackletTree-'$1'TeV-'$3'-2.root", '${_MC_SPLIT[2]}', '${_MC_SPLIT[3]}')' &> logs/${1}tev_mc_2.log && root -l -b -q 'analyze_trackletTree.C+("'$_MC_PTREE'", "TrackletTree-'$1'TeV-'$3'-3.root", '${_MC_SPLIT[3]}', '${_MC_SPLIT[4]}')' &> logs/${1}tev_mc_3.log &

    root -l -b -q 'analyze_trackletTree.C+("'$_MC_PTREE'", "TrackletTree-'$1'TeV-'$3'-4.root", '${_MC_SPLIT[4]}', '${_MC_SPLIT[5]}')' &> logs/${1}tev_mc_4.log && root -l -b -q 'analyze_trackletTree.C+("'$_MC_PTREE'", "TrackletTree-'$1'TeV-'$3'-5.root", '${_MC_SPLIT[5]}', '${_MC_SPLIT[6]}')' &> logs/${1}tev_mc_5.log &
    wait
fi

if [ $2 -ne 2 ]; then
    ./configure.sh $1 0

    echo "producing data tracklettree..."
    root -l -b -q 'analyze_trackletTree.C++("'$_DATA_PTREE'", "TrackletTree-COMPILE-DUMMY.root", 0, 1)'
    rm TrackletTree-COMPILE-DUMMY.root

    root -l -b -q 'analyze_trackletTree.C+("'$_DATA_PTREE'", "TrackletTree-'$_RUN'-'$3'-0.root", '${_DATA_SPLIT[0]}', '${_DATA_SPLIT[1]}')' &> logs/${1}tev_data_0.log && root -l -b -q 'analyze_trackletTree.C+("'$_DATA_PTREE'", "TrackletTree-'$_RUN'-'$3'-1.root", '${_DATA_SPLIT[1]}', '${_DATA_SPLIT[2]}')' &> logs/${1}tev_data_1.log &

    root -l -b -q 'analyze_trackletTree.C+("'$_DATA_PTREE'", "TrackletTree-'$_RUN'-'$3'-2.root", '${_DATA_SPLIT[2]}', '${_DATA_SPLIT[3]}')' &> logs/${1}tev_data_2.log && root -l -b -q 'analyze_trackletTree.C+("'$_DATA_PTREE'", "TrackletTree-'$_RUN'-'$3'-3.root", '${_DATA_SPLIT[3]}', '${_DATA_SPLIT[4]}')' &> logs/${1}tev_data_3.log &

    root -l -b -q 'analyze_trackletTree.C+("'$_DATA_PTREE'", "TrackletTree-'$_RUN'-'$3'-4.root", '${_DATA_SPLIT[4]}', '${_DATA_SPLIT[5]}')' &> logs/${1}tev_data_4.log && root -l -b -q 'analyze_trackletTree.C+("'$_DATA_PTREE'", "TrackletTree-'$_RUN'-'$3'-5.root", '${_DATA_SPLIT[5]}', '${_DATA_SPLIT[6]}')' &> logs/${1}tev_data_5.log &
    wait
fi
