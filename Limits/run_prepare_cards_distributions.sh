#!/bin/bash -x

#TYPE="llstau_maximally-mixed"
#TYPE="llstau_mass-degenerate"

SUFFIX="$1"
CHANNEL="$2"
ARGS="$3"

if [ -n "$SUFFIX" ]; then
    SUFFIX="_${SUFFIX}"
fi

python3 prepare_cards.py --configs \
configs/distributions/2016_preVFP/config_distribution_${CHANNEL}.yaml \
configs/distributions/2016_postVFP/config_distribution_${CHANNEL}.yaml \
configs/distributions/2017/config_distribution_${CHANNEL}.yaml \
configs/distributions/2018/config_distribution_${CHANNEL}.yaml \
--yields_uncs \
--yields_uncs_sigs configs/distributions/sig_list_for-distributions.txt \
--nocards \
--outdir results/distributions${SUFFIX} \
${ARGS}
