#!/bin/bash

./run-final-result.sh 5 latest
./run-final-result.sh 8 latest

./run-systematics.sh 5
./run-systematics.sh 8

g++ calc_systematics.C -o calc_systematics -O2 -Wall `root-config --cflags --libs`
./calc_systematics rootfiles/merged-285090-latest.root 5tev.list 5tev
./calc_systematics rootfiles/merged-285832-latest.root 8tev.list 8tev

./make_final_plots results.list rootfiles/results.root gen.list
