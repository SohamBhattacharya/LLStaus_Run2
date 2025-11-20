#!/bin/bash -x

# Exit if a command fails
set -Eeu -o pipefail

SUFFIX="$1"
ARGS="$2"

if [ -n "$SUFFIX" ]; then
    SUFFIX="_${SUFFIX}"
fi

DISTRIBUTIONS=(
    "jet2_pt"
    "METpt"
    "mt2_j1_j2_MET"
)

REGIONS=(
    #"BRT1"
    "BRT2"
)

ERAS=(
    #"2016_preVFP"
    #"2016_postVFP"
    #"2017"
    #"2018"
    "2016_preVFP.2016_postVFP.2017.2018"
)

declare -A LUMITEXTS=(
    ["2016_preVFP"]="19.5 fb^{-1} (13 TeV)"
    ["2016_postVFP"]="16.8 fb^{-1} (13 TeV)"
    ["2017"]="41.5 fb^{-1} (13 TeV)"
    ["2018"]="59.8 fb^{-1} (13 TeV)"
    ["2016_preVFP.2016_postVFP.2017.2018"]="138 fb^{-1} (13 TeV)"
)

LIMITSDIR=distributions${SUFFIX}
INDIR=results/${LIMITSDIR}/llstau_maximally-mixed/yields_and_systematics
OUTDIR=${INDIR}/figures_and_tables

for dist in ${DISTRIBUTIONS[@]}; do
    for region in ${REGIONS[@]}; do
        for era in ${ERAS[@]}; do
            ./yaml_to_table.py \
            --type yields \
            --plot \
            --lumitext "${LUMITEXTS[${era}]}" \
            --config configs/tables/distributions/config_distribution_${dist}_${region}.yaml \
            --input ${INDIR}/yields_channels_${region}_${dist}_eras_${era}.yaml:${era} \
            --outdir $OUTDIR \
            ${ARGS}
            
            #./yaml_to_table.py \
            #--type systematics \
            #--config configs/tables/distributions/config_systematics_table_${region}.yaml \
            #--input ${INDIR}/systematics_channels_${region}_eras_${era}.yaml:${era} \
            #--outdir $OUTDIR \
            #${ARGS}
            
            #printf "\n\n"
        done
    done
done
