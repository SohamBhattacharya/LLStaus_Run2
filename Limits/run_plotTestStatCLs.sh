#!/bin/bash

DIR="$1"

pushd $DIR

$CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/plotTestStatCLs.py \
--input higgsCombine.Test.HybridNew.mH120.merged.root \
--poi r \
--val all \
--mass 120 \
--save-as-pdf \
-E \
-q 0.975

popd
