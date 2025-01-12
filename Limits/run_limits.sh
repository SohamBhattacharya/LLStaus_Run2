#!/bin/bash

#set -f
#PATTERN=$1
#echo $PATTERN
#set +f
#combineTool.py -M T2W -i "$PATTERN" -o workspace.root --parallel 15
#combineTool.py -M Asymptotic -t -1 -d "$PATTERN"/workspace.root --there -n .limit --parallel 15

DIR=results/limits/llstau_maximally-mixed/channels_all/eras_all
#DIR=results/limits/llstau_mass-degenerate/channels_all/eras_all

nice -n 10 combineTool.py -M T2W -i $DIR/SMS-TStauStau_MStau-*/card_*.txt -o workspace.root --parallel 20

#nice -n 10 combineTool.py -M AsymptoticLimits -t -1 -d $DIR/SMS-TStauStau_MStau-*/workspace.root --there -n .limit --parallel 20
nice -n 10 combineTool.py -M AsymptoticLimits -d $DIR/SMS-TStauStau_MStau-*/workspace.root --there -n .limit --parallel 20

nice -n 10 combineTool.py -M CollectLimits $DIR/SMS-TStauStau_MStau-*/*.limit.* --use-dirs -o $DIR/limits.json
