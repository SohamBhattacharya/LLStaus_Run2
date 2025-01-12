#!/bin/bash

DIR=$1

pushd $1

FLIST=( \
GenStau1_pathL.pdf GenTau1_vertexR.pdf Jet1_dxy.pdf Jet1_disTauTag_score1.pdf \
GenStau2_pathL.pdf GenTau2_vertexR.pdf Jet2_dxy.pdf Jet2_disTauTag_score1.pdf \
)

echo "${FLIST[@]}"

pdfjam --suffix nup --delta '-0.3cm 0cm' --landscape --nup 4x2 --outfile merge.pdf "${FLIST[@]}"
pdfcrop merge.pdf merge.pdf
pdftoppm -png -r 600 -cropbox merge.pdf merge

popd

echo "Done"
