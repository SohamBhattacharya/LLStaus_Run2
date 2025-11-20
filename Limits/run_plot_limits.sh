#!/bin/bash -x

set -eEu

TYPE="$1"
SUFFIX="$2"
ARGS="$3"

if [ -n "$SUFFIX" ]; then
    SUFFIX="_${SUFFIX}"
fi

if [ "$TYPE" = "mm" ]; then
    TYPE="maximally-mixed"
    TEXT="Maximally mixed scenario"
elif [ "$TYPE" = "md" ]; then
    TYPE="mass-degenerate"
    TEXT="Mass degenerate scenario"
else
    echo "Invalid TYPE ${TYPE}; must be mm or md"
    exit 1
fi


#DIR=~/mnt/desy_dust/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/results/limits${SUFFIX}/llstau_${TYPE}/channels_all/eras_all
#DIR=~/mnt/desy_dust/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/results/limits_test-distau-syst-20percent/llstau_${TYPE}/channels_all/eras_all
#DIR=~/mnt/desy_dust/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/results/limits_dxy-gt-0p2_8-bins/llstau_${TYPE}/channels_all/eras_all
#DIR=results/limits${SUFFIX}/llstau_${TYPE}/channels_all/eras_all
DIR=results/limits${SUFFIX}/llstau_${TYPE}/channels_BRT2/eras_all

XSECFILE=../Analysis/configs/crosssections_stau_${TYPE}_hepi-fast.csv

./plot_limits.py \
--jsons $DIR/SMS-TStauStau_MStau-*/limits/limits.json \
--xsecfile $XSECFILE \
--extratext "$TEXT" \
--cmsextratext "Preliminary" \
--outdir $DIR/limits \
${ARGS}
