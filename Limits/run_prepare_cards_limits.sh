#!/bin/bash -x

TYPE="llstau_maximally-mixed"
#TYPE="llstau_mass-degenerate"

ARGS="$1"
OUTSUFFIX="$2"

python3 prepare_cards.py --configs \
configs/limits/$TYPE/2016_preVFP/config_datacard_BRT2.yaml \
configs/limits/$TYPE/2016_postVFP/config_datacard_BRT2.yaml \
configs/limits/$TYPE/2017/config_datacard_BRT2.yaml \
configs/limits/$TYPE/2018/config_datacard_BRT2.yaml \
--yields_uncs \
--yields_uncs_sigs configs/limits/sig_list_for-tables.txt \
--outdir results/limits${OUTSUFFIX} \
--combpars era \
$1
