#!/bin/bash

# Examples
# ./get-xsecs_hepi-fast.sh ../Analysis/configs/13000_sleptons_1000015_-1000015_NNLL.json ../Analysis/configs/crosssections_stau_maximally-mixed_hepi-fast.csv
# ./get-xsecs_hepi-fast.sh ../Analysis/configs/pp13_stau_LR_NLO+NLL_PDF4LHC.json ../Analysis/configs/crosssections_stau_mass-degenerate_hepi-fast.csv

JSON=$1
OUTPUT=$2
MASSLIST="configs/limits/sig_mass_list.txt"

if [ -z "${JSON}" ] || [ -z "${OUTPUT}" ]; then
    echo "Error. Usage:"
    echo "get-xsecs_hepi-fast.sh [INPUT JSON] [OUTPUT]"
    exit 1
fi

#hepi-fast --help
OUTPUT_TMP="/tmp/tmp_hepi-fast.txt"
HEADER="#Mass [GeV] , Central value [fb] , error down [fb] , error up [fb] , error pdf down [fb] , error pdf up [fb] , error scale down [fb] , error scale up [fb]"

hepi-fast $JSON < $MASSLIST > $OUTPUT_TMP

echo $HEADER > $OUTPUT

# Exclude the 1st column -- it is just 0
# Convert pb to fb (*1000)

awk '{ print \
    $2 " , " \
    $3*1000 " , " \
    $4*1000 " , " \
    $5*1000 " , " \
    $6*1000 " , " \
    $7*1000 " , " \
    $8*1000 " , " \
    $9*1000       \
}' $OUTPUT_TMP >> $OUTPUT