#!/bin/bash

MASSES=(
    "MStau-100_ctau-50mm"
    "MStau-100_ctau-100mm"
    "MStau-200_ctau-50mm"
    "MStau-200_ctau-100mm"
)

declare -A LIMTYPES=(
    ["obs"]=""
    ["exp2p5"]="--expectedFromGrid=0.025"
    ["exp16"]="--expectedFromGrid=0.16"
    ["exp50"]="--expectedFromGrid=0.5"
    ["exp84"]="--expectedFromGrid=0.84"
    ["exp97p5"]="--expectedFromGrid=0.975"
)

DIR="results/limits_signal_v12/llstau_maximally-mixed/channels_BRT2/eras_all"

for mass in ${MASSES[@]}; do
    for lim in ${LIMTYPES[@]}; do
        lim_suffix=""
        if [ -n "$lim" ]; then
            lim_suffix="_${lim}"
        fi
        
        wspace="${DIR}/SMS-TStauStau_MStau-200_ctau-100mm_mLSP-1/limits"
        
    done