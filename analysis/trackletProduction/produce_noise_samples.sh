#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: ./produce_noise_samples [5/8]"
    exit 1
fi

case $1 in
    5)
        echo "producing monte carlo noise samples for 5.02 TeV"
        _MC_PTREE="/data/biran/pixeltrees/PixelTree-EPOS-5TeV-HLT.root"
        _MC_SPLIT=(0 80000 160000 240000 320000 400000 480000)
        ;;
    8)
        echo "producing monte carlo noise samples for 8.16 TeV"
        _MC_PTREE="/data/biran/pixeltrees/PixelTree-EPOS-8TeV-HLT.root"
        _MC_SPLIT=(0 135000 270000 405000 540000 675000 810000)
        ;;
    *)
        echo "check sqrt(snn) energy!"
        exit 1
        ;;
esac

./configure.sh $1 1

echo "compilation test..."
root -l -b -q 'analyze_trackletTree.C++("'$_MC_PTREE'", "TrackletTree-COMPILE-DUMMY.root", 0, 1)'
rm TrackletTree-COMPILE-DUMMY.root

echo "producing sample with additional noise..."
root -l -b -q 'analyze_trackletTree.C+("'$_MC_PTREE'", "TrackletTree-MCSYS-'${1}'-NOISE-0.root", '${_MC_SPLIT[0]}', '${_MC_SPLIT[1]}', 0, 1, 0, 1, 1, 0.05, 0.05, 0.05)' &> logs/${1}tev_noise_0_mc.log && root -l -b -q 'analyze_trackletTree.C+("'$_MC_PTREE'", "TrackletTree-MCSYS-'${1}'-NOISE-1.root", '${_MC_SPLIT[1]}', '${_MC_SPLIT[2]}', 0, 1, 0, 1, 1, 0.05, 0.05, 0.05)' &> logs/${1}tev_noise_1_mc.log &

root -l -b -q 'analyze_trackletTree.C+("'$_MC_PTREE'", "TrackletTree-MCSYS-'${1}'-NOISE-2.root", '${_MC_SPLIT[2]}', '${_MC_SPLIT[3]}', 0, 1, 0, 1, 1, 0.05, 0.05, 0.05)' &> logs/${1}tev_noise_2_mc.log && root -l -b -q 'analyze_trackletTree.C+("'$_MC_PTREE'", "TrackletTree-MCSYS-'${1}'-NOISE-3.root", '${_MC_SPLIT[3]}', '${_MC_SPLIT[4]}', 0, 1, 0, 1, 1, 0.05, 0.05, 0.05)' &> logs/${1}tev_noise_3_mc.log &

root -l -b -q 'analyze_trackletTree.C+("'$_MC_PTREE'", "TrackletTree-MCSYS-'${1}'-NOISE-4.root", '${_MC_SPLIT[4]}', '${_MC_SPLIT[5]}', 0, 1, 0, 1, 1, 0.05, 0.05, 0.05)' &> logs/${1}tev_noise_4_mc.log && root -l -b -q 'analyze_trackletTree.C+("'$_MC_PTREE'", "TrackletTree-MCSYS-'${1}'-NOISE-5.root", '${_MC_SPLIT[5]}', '${_MC_SPLIT[6]}', 0, 1, 0, 1, 1, 0.05, 0.05, 0.05)' &> logs/${1}tev_noise_5_mc.log &
wait
