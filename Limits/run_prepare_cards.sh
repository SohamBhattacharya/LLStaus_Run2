#!/bin/bash

TYPE="llstau_maximally-mixed"
#TYPE="llstau_mass-degenerate"

python3 prepare_cards.py --configs \
configs/limits/$TYPE/2016_preVFP/config_datacard_BRT2.yaml \
configs/limits/$TYPE/2016_postVFP/config_datacard_BRT2.yaml \
configs/limits/$TYPE/2017/config_datacard_BRT2.yaml \
configs/limits/$TYPE/2018/config_datacard_BRT2.yaml \
--outdir results/limits
