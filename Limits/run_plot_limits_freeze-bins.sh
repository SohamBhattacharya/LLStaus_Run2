#!/bin/bash -x

SUFFIX="$1"
MODE="$2"

if [ -n "$SUFFIX" ]; then
    SUFFIX="_${SUFFIX}"
fi

TYPE="maximally-mixed"
TEXT="Maximally mixed scenario"

#TYPE="mass-degenerate"
#TEXT="Mass degenerate scenario"

DIR=results/limits${SUFFIX}/llstau_${TYPE}/channels_BRT2/eras_all

XSECFILE=../Analysis/configs/crosssections_stau_${TYPE}_hepi-fast.csv

BINS=(
    "bin3_BRT2"
    "bin4_BRT2"
    "bin5_BRT2"
    "bin6_BRT2"
    "bin7_BRT2"
    "bin8_BRT2"
    "bin9_BRT2"
    "bin10_BRT2"
)

for bin in ${BINS[@]}; do
    
    bin_suffix="_${MODE}_${bin}"
    
    ./plot_limits.py \
    --jsons ${DIR}/SMS-TStauStau_MStau-*/limits${bin_suffix}/limits.json \
    --xsecfile $XSECFILE \
    --extratext "$TEXT" \
    --cmsextratext "Private Work" \
    --output $DIR/limits${bin_suffix}/limits${bin_suffix}.root
done
