#!/usr/bin/env bash

N=$1
Q=$2
SD=$3
M=$4
BETA=$5
TOURS=$6
ENUMGHBOUND=$7
PRNGSEED=$8
SE93=$9

# # generating lattices using fplll
# ./fplll/out/bin/latticegen -randseed 1234 u 20 4 | sage bkz.py -b 10 -tree-filename tree.json -tree-stats-filename stats.json
# ./fplll/out/bin/latticegen -randseed 1234 u 20 4 | /home/fernando/source/QEnumeration/Code/EnumerationCandidates/fplll/out/bin/fplll -a bkz -v -b 10 -bkzautoabort -bkzmaxloops 20 -bkzghbound 1.1 -bkzdumpenumtrees tree.json -bkzdumpenumtreestats stats.json
# sage plot_children_predictions.py -tree-stats-filename stats.json


# generate LWE lattices

if [[ -z "$N" ]]; then N=72; fi
if [[ -z "$Q" ]]; then Q=97; fi
if [[ -z "$SD" ]]; then SD=1.0; fi
if [[ -z "$M" ]]; then M=87; fi
if [[ -z "$BETA" ]]; then BETA=30; fi
if [[ -z "$TOURS" ]]; then TOURS=20; fi
if [[ -z "$ENUMGHBOUND" ]]; then ENUMGHBOUND=1.1; fi
if [[ -z "$PRNGSEED" ]]; then PRNGSEED=1337; fi

STATSFN=./stats/stats-n$N-q$Q-sd$SD-m$M-beta$BETA-tours$TOURS-ghbound$ENUMGHBOUND-seed$PRNGSEED.json
PLOTDIR=./plots/n$N-q$Q-sd$SD-m$M-beta$BETA-tours$TOURS-ghbound$ENUMGHBOUND-seed$PRNGSEED/

# echo N = $N
# echo Q = $Q
# echo SD = $SD
# echo M = $M
# echo BETA = $BETA
# echo TOURS = $TOURS
# echo ENUMGHBOUND = $ENUMGHBOUND
# echo PRNGSEED = $PRNGSEED
# echo STATSFN = $STATSFN
# exit 1

if [[ -z "$SE93" ]]
then
    # STRAT=fplll/strategies/default.json
    STRAT=no-pruning.json
    sage lwe.py -randseed ${PRNGSEED} -n ${N} -q ${Q} -sd ${SD} -m ${M} | ./out/bin/fplll -a bkz -s ${STRAT}  -v -b ${BETA} -bkzautoabort -bkzmaxloops ${TOURS} -bkzghbound ${ENUMGHBOUND} -bkzdumpenumtreestats ${STATSFN} || true
    # get out secret info in LWE lattice
    # sage lwe.py -randseed 1337 -n 72 -q 97 -sd 1. -m 87 -printsol
    sage plot_children_predictions.py -n ${N} -q ${Q} -sd ${SD} -m ${M} -beta ${BETA} -tours ${TOURS} -ghfactor ${ENUMGHBOUND} -tree-stats-filename ${STATSFN} -plot-dir ${PLOTDIR}
else
    sage lwe.py -randseed ${PRNGSEED} -n ${N} -q ${Q} -sd ${SD} -m ${M} | sage bkz.py -SE93 -fpylll -verbose -b ${BETA} -maxtours ${TOURS}
fi


