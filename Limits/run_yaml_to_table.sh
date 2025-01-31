#!/bin/bash -x

# Exit if a command fails
set -Eeu -o pipefail

ERAS=(
    "2016_preVFP"
    "2016_postVFP"
    "2017"
    "2018"
    "2016_preVFP.2016_postVFP.2017.2018"
)

INDIR=/home/soham/mnt/desy_dust/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/results/limits_blinded_signal_v6/llstau_maximally-mixed
OUTDIR=results/limits_blinded_signal_v6/llstau_maximally-mixed

for era in ${ERAS[@]}; do
    ./yaml_to_table.py --type yields --config configs/tables/limits/config_yields_table.yaml \
    --input ${INDIR}/yields_channels_BRT2_eras_${era}.yaml:${era} \
    --outdir $OUTDIR
    
    ./yaml_to_table.py --type systematics --config configs/tables/limits/config_systematics_table.yaml \
    --input ${INDIR}/systematics_channels_BRT2_eras_${era}.yaml:${era} \
    --outdir $OUTDIR
    
    #printf "\n\n"
done
