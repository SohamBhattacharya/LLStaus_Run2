#!/bin/bash

set -Eeu -o pipefail
#set -x

INDIR="$1"

CTAUS=(
    1mm
    2mm
    3mm
    4mm
    5mm
    6p5mm
    8p5mm
    10mm
    20mm
    30mm
    40mm
    50mm
    60mm
    80mm
    100mm
    150mm
    200mm
    300mm
    400mm
    500mm
    700mm
)

for ctau in "${CTAUS[@]}"; do
    outfile="${INDIR}/compare_limits_${ctau}_merged.png"
    outfile_basename=$(basename "$outfile")
    
    montage -tile 3x3 -mode concatenate $(find $INDIR/*_${ctau}.png | grep -e obs -e exp | grep -v $outfile_basename| sort -V) $outfile
    echo "Created combined file for ${ctau}: ${outfile}"
done