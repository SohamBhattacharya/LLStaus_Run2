
#!/bin/bash

# Exit if a command fails
set -Eeu -o pipefail

# Prints command before executing
set -o xtrace

ERA="$1"
SUFFIX="$2"

./run_prepare_cards_DisTauSF_mass.py --era $ERA --outsuffix "$SUFFIX"

if [ -n "$SUFFIX" ]; then
    SUFFIX="_${SUFFIX}"
fi

./fit_DisTauSF.py \
--wspace \
--fit \
--collfits \
--postfit \
--scan \
--impacts \
--indir "results/DisTauSF/DisTauSF_mass_nbins1-60-80${SUFFIX}/DisTauSF/channels_all/eras_${ERA}/ZMT_wp-p*"

