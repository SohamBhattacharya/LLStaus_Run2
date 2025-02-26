#!/bin/bash -x

# Exit if a command fails
set -Eeu -o pipefail

SUFFIX="$1"
ARGS="$2"

if [ -n "$SUFFIX" ]; then
    SUFFIX="_${SUFFIX}"
fi

REGIONS=(
#    "BRT1"
    "BRT2"
)

ERAS=(
    "2016_preVFP"
    "2016_postVFP"
    "2017"
    "2018"
    "2016_preVFP.2016_postVFP.2017.2018"
)

declare -A LUMITEXTS=(
    ["2016_preVFP"]="19.5 fb^{-1} (13 TeV)"
    ["2016_postVFP"]="16.8 fb^{-1} (13 TeV)"
    ["2017"]="41.5 fb^{-1} (13 TeV)"
    ["2018"]="59.8 fb^{-1} (13 TeV)"
    ["2016_preVFP.2016_postVFP.2017.2018"]="138 fb^{-1} (13 TeV)"
)

#LIMITSDIR=limits_signal_v9_blinded
LIMITSDIR=limits${SUFFIX}
#INDIR=/home/soham/mnt/desy_dust/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/results/${LIMITSDIR}/llstau_maximally-mixed
INDIR=results/${LIMITSDIR}/llstau_maximally-mixed/yields_and_systematics
#OUTDIR=results/${LIMITSDIR}/llstau_maximally-mixed/yields_and_systematics
OUTDIR=${INDIR}/figures_and_tables

for region in ${REGIONS[@]}; do
    for era in ${ERAS[@]}; do
        ./yaml_to_table.py \
        --type yields \
        --plot \
        --lumitext "${LUMITEXTS[${era}]}" \
        --config configs/tables/limits/config_yields_table_${region}.yaml \
        --input ${INDIR}/yields_channels_${region}_eras_${era}.yaml:${era} \
        --outdir $OUTDIR \
        ${ARGS}
        
        ./yaml_to_table.py \
        --type systematics \
        --config configs/tables/limits/config_systematics_table_${region}.yaml \
        --input ${INDIR}/systematics_channels_${region}_eras_${era}.yaml:${era} \
        --outdir $OUTDIR \
        ${ARGS}
        
        #printf "\n\n"
    done
done
