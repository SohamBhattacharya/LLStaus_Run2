#!/bin/bash -x

set -eEu

SUFFIX="$1"
ARGS="$2"

if [ -n "$SUFFIX" ]; then
    SUFFIX="_${SUFFIX}"
fi

python3 prepare_cards.py --configs \
configs/compare-predictions/2016_preVFP/config_datacard_BRT2_compare-predictions.yaml \
configs/compare-predictions/2016_postVFP/config_datacard_BRT2_compare-predictions.yaml \
configs/compare-predictions/2017/config_datacard_BRT2_compare-predictions.yaml \
configs/compare-predictions/2018/config_datacard_BRT2_compare-predictions.yaml \
--yields_uncs \
--outdir results/compare-predictions${SUFFIX} \
--combpars era \
--nocards \
${ARGS}
