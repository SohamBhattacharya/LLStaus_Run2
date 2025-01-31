#!/bin/bash -x

SUFFIX="$1"

TYPE="maximally-mixed"
TEXT="Maximally mixed scenario"

#TYPE="mass-degenerate"
#TEXT="Mass degenerate scenario"

#DIR=~/mnt/desy_dust/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/results/limits${SUFFIX}/llstau_${TYPE}/channels_all/eras_all
#DIR=~/mnt/desy_dust/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/results/limits_test-distau-syst-20percent/llstau_${TYPE}/channels_all/eras_all
#DIR=~/mnt/desy_dust/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/results/limits_dxy-gt-0p2_8-bins/llstau_${TYPE}/channels_all/eras_all
DIR=results/limits${SUFFIX}/llstau_${TYPE}/channels_all/eras_all

XSECFILE=../Analysis/configs/crosssections_stau_${TYPE}_hepi-fast.csv

./plot_limits.py --jsons $DIR/SMS-TStauStau_MStau-*/limits/limits.json --xsecfile $XSECFILE --extratext "$TEXT" --cmsextratext "Private Work" --output $DIR/limits.root
