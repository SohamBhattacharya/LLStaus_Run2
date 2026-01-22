#!/bin/bash -x

# Exit if a command fails
set -Eeu -o pipefail

INDIR="$1"

#SUFFIX="$2"
#
#if [ -n "$SUFFIX" ]; then
#    SUFFIX="_${SUFFIX}"
#fi

FITTYPES=(
    #"s"
    "b"
)

MASSES=(
    "MStau-100_ctau-50mm"
    #"MStau-100_ctau-100mm"
    #"MStau-200_ctau-50mm"
    #"MStau-200_ctau-100mm"
)

ERAS=(
    #"2016_preVFP"
    #"2016_postVFP"
    #"2017"
    #"2018"
    
    "added"
)

OUTDIR="${INDIR}/prefit-postfit"

for mass in ${MASSES[@]}; do
    for fit in ${FITTYPES[@]}; do
        
        dir="${INDIR}/SMS-TStauStau_${mass}_mLSP-1/fit-diagnostics"
        
        ./plot_prefit-postfit_limits.py \
        --input ${dir}/postfit_${fit}.root \
        --era ${ERAS[@]} \
        --channels BRT2 \
        --outdir "${OUTDIR}/${fit}" \
        --type ${fit}
        
        
        ./plot_prefit-postfit_limits.py \
        --input ${dir}/postfit_${fit}.root \
        --era ${ERAS[@]} \
        --channels BRT2 \
        --outdir "${OUTDIR}/${fit}_nosig" \
        --type ${fit} \
        --nosig \
        --cmsextratext ""

    done
done