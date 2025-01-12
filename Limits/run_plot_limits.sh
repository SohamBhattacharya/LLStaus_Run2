#!/bin/bash

DIR=~/mnt/nfs_dust/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/results/limits/llstau_maximally-mixed/channels_all/eras_all
#DIR=~/mnt/nfs_dust/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/results/limits/llstau_mass-degenerate/channels_all/eras_all

./plot_limits.py --jsons $DIR/*.json --output $DIR/limits.root
