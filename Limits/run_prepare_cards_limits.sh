#!/bin/bash -x

set -eEu

#TYPE="llstau_maximally-mixed"
TYPE="llstau_mass-degenerate"

ARGS="$1"
SUFFIX="$2"
CHANNEL="$3"

if [ -n "$SUFFIX" ]; then
    SUFFIX="_${SUFFIX}"
fi

python3 prepare_cards.py --configs \
configs/limits/$TYPE/2016_preVFP/config_datacard_${CHANNEL}.yaml \
configs/limits/$TYPE/2016_postVFP/config_datacard_${CHANNEL}.yaml \
configs/limits/$TYPE/2017/config_datacard_${CHANNEL}.yaml \
configs/limits/$TYPE/2018/config_datacard_${CHANNEL}.yaml \
--yields_uncs \
--yields_uncs_sigs configs/limits/sig_list_for-tables.txt \
--outdir results/limits${SUFFIX} \
--combpars era \
${ARGS}
