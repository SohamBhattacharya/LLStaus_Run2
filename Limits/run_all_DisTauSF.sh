
#!/bin/bash

# Exit if a command fails
set -Eeu -o pipefail

# Prints command before executing
set -o xtrace

ERA="$1"

./run_prepare_cards_DisTauSF_mass.py --era $ERA

./fit_DisTauSF.py \
--wspace \
--fit \
--collfits \
--postfit \
--scan \
--impacts \
--indir "results/DisTauSF/DisTauSF_mass_nbins1-60-80/DisTauSF/channels_all/eras_${ERA}/ZMT_wp-p*"

