#!/bin/bash -x

# Exit if a command fails
set -Eeu -o pipefail

SUFFIX="$1"
ARGS="$2"

if [ -n "$SUFFIX" ]; then
    SUFFIX="_${SUFFIX}"
fi

REGIONS=(
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

INDIR=results/compare-predictions${SUFFIX}/llstau_maximally-mixed/yields_and_systematics
OUTDIR=${INDIR}/figures_and_tables

for era in ${ERAS[@]}; do
    ./yaml_to_table.py \
    --type yields \
    --plot \
    --lumitext "${LUMITEXTS[${era}]}" \
    --config configs/tables/compare-predictions/config_yields_table_BRT2_compare-predictions.yaml \
    --input ${INDIR}/yields_channels_BRT2_eras_${era}.yaml:${era} \
    --outdir $OUTDIR \
    ${ARGS}
done
