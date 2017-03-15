#!/bin/bash

if [ $# -lt 2 ]; then
    echo "Usage: ./produce_accep_samples [5/8] [0/1/2]"
    echo "       0: ALL, 1: DATA, 2: MC"
    exit 1
fi

case $1 in
    5)
        echo "producing acceptance samples for 5.02 TeV"

        _RUN=285090
        _MC_PTREE="/data/biran/pixeltrees/PixelTree-EPOS-5TeV-HLT.root"
        _DATA_PTREE="/data/biran/pixeltrees/2016pPb/Express/PixelTree-285090.root"
        _MC_SPLIT=(0 250000 500000 0 250000 500000 0 250000 500000 0 250000 500000)
        _DATA_SPLIT=(0 960000 1920000 0 960000 1920000 0 960000 1920000 0 960000 1920000)
        ;;
    8)
        echo "producing acceptance samples for 8.16 TeV"

        _RUN=285832
        _MC_PTREE="/data/biran/pixeltrees/PixelTree-EPOS-8TeV-HLT.root"
        _DATA_PTREE="/data/biran/pixeltrees/2016pPb/Express/PixelTree-285832.root"
        _MC_SPLIT=(0 202000 404000 404000 606000 808000 0 202000 404000 404000 606000 808000)
        _DATA_SPLIT=(3000000 3375000 3750000 3750000 4125000 4500000 4500000 4875000 5250000 5250000 5625000 6000000)
        ;;
    *)
        echo "check sqrt(snn) energy!"
        exit 1
        ;;
esac

if [ $2 -ne 1 ]; then
    echo "producing mc acceptance samples"
    ./configure.sh $1 1

    root -l -b -q 'analyze_trackletTree.C++("'$_MC_PTREE'", "TrackletTree-COMPILE-DUMMY.root", 0, 1)'
    rm TrackletTree-COMPILE-DUMMY.root

    root -l -b -q 'analyze_trackletTree.C+("'$_MC_PTREE'", "TrackletTree-'$1'TeV-ACCEP-0.root", '${_MC_SPLIT[0]}', '${_MC_SPLIT[1]}', 1, 1, 0, 0, 0)' &> logs/${1}tev_mc_accep_0.log && root -l -b -q 'analyze_trackletTree.C+("'$_MC_PTREE'", "TrackletTree-'$1'TeV-ACCEP-1.root", '${_MC_SPLIT[1]}', '${_MC_SPLIT[2]}', 1, 1, 0, 0, 0)' &> logs/${1}tev_mc_accep_1.log &

    root -l -b -q 'analyze_trackletTree.C+("'$_MC_PTREE'", "TrackletTree-'$1'TeV-ACCEP-2.root", '${_MC_SPLIT[3]}', '${_MC_SPLIT[4]}', 1, 1, 0, 0, 0)' &> logs/${1}tev_mc_accep_2.log && root -l -b -q 'analyze_trackletTree.C+("'$_MC_PTREE'", "TrackletTree-'$1'TeV-ACCEP-3.root", '${_MC_SPLIT[4]}', '${_MC_SPLIT[5]}', 1, 1, 0, 0, 0)' &> logs/${1}tev_mc_accep_3.log &

    root -l -b -q 'analyze_trackletTree.C+("'$_MC_PTREE'", "TrackletTree-'$1'TeV-ACCEP-4.root", '${_MC_SPLIT[6]}', '${_MC_SPLIT[7]}', 1, 1, 0, 0, 0)' &> logs/${1}tev_mc_accep_4.log && root -l -b -q 'analyze_trackletTree.C+("'$_MC_PTREE'", "TrackletTree-'$1'TeV-ACCEP-5.root", '${_MC_SPLIT[7]}', '${_MC_SPLIT[8]}', 1, 1, 0, 0, 0)' &> logs/${1}tev_mc_accep_5.log &

    root -l -b -q 'analyze_trackletTree.C+("'$_MC_PTREE'", "TrackletTree-'$1'TeV-ACCEP-6.root", '${_MC_SPLIT[9]}', '${_MC_SPLIT[10]}', 1, 1, 0, 0, 0)' &> logs/${1}tev_mc_accep_6.log && root -l -b -q 'analyze_trackletTree.C+("'$_MC_PTREE'", "TrackletTree-'$1'TeV-ACCEP-7.root", '${_MC_SPLIT[10]}', '${_MC_SPLIT[11]}', 1, 1, 0, 0, 0)' &> logs/${1}tev_mc_accep_7.log &
    wait
fi

if [ $2 -ne 2 ]; then
    echo "producing data acceptance samples..."
    ./configure.sh $1 0

    root -l -b -q 'analyze_trackletTree.C++("'$_DATA_PTREE'", "TrackletTree-COMPILE-DUMMY.root", 0, 1)'
    rm TrackletTree-COMPILE-DUMMY.root

    root -l -b -q 'analyze_trackletTree.C+("'$_DATA_PTREE'", "TrackletTree-'$_RUN'-ACCEP-0.root", '${_DATA_SPLIT[0]}', '${_DATA_SPLIT[1]}', 1, 1, 0, 0, 0)' &> logs/${1}tev_data_accep_0.log && root -l -b -q 'analyze_trackletTree.C+("'$_DATA_PTREE'", "TrackletTree-'$_RUN'-ACCEP-1.root", '${_DATA_SPLIT[1]}', '${_DATA_SPLIT[2]}', 1, 1, 0, 0, 0)' &> logs/${1}tev_data_accep_1.log &

    root -l -b -q 'analyze_trackletTree.C+("'$_DATA_PTREE'", "TrackletTree-'$_RUN'-ACCEP-2.root", '${_DATA_SPLIT[3]}', '${_DATA_SPLIT[4]}', 1, 1, 0, 0, 0)' &> logs/${1}tev_data_accep_2.log && root -l -b -q 'analyze_trackletTree.C+("'$_DATA_PTREE'", "TrackletTree-'$_RUN'-ACCEP-3.root", '${_DATA_SPLIT[4]}', '${_DATA_SPLIT[5]}', 1, 1, 0, 0, 0)' &> logs/${1}tev_data_accep_3.log &

    root -l -b -q 'analyze_trackletTree.C+("'$_DATA_PTREE'", "TrackletTree-'$_RUN'-ACCEP-4.root", '${_DATA_SPLIT[6]}', '${_DATA_SPLIT[7]}', 1, 1, 0, 0, 0)' &> logs/${1}tev_data_accep_4.log && root -l -b -q 'analyze_trackletTree.C+("'$_DATA_PTREE'", "TrackletTree-'$_RUN'-ACCEP-5.root", '${_DATA_SPLIT[7]}', '${_DATA_SPLIT[8]}', 1, 1, 0, 0, 0)' &> logs/${1}tev_data_accep_5.log &

    root -l -b -q 'analyze_trackletTree.C+("'$_DATA_PTREE'", "TrackletTree-'$_RUN'-ACCEP-6.root", '${_DATA_SPLIT[9]}', '${_DATA_SPLIT[10]}', 1, 1, 0, 0, 0)' &> logs/${1}tev_data_accep_6.log && root -l -b -q 'analyze_trackletTree.C+("'$_DATA_PTREE'", "TrackletTree-'$_RUN'-ACCEP-7.root", '${_DATA_SPLIT[10]}', '${_DATA_SPLIT[11]}', 1, 1, 0, 0, 0)' &> logs/${1}tev_data_accep_7.log &
    wait
fi
