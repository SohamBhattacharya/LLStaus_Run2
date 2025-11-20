#!/bin/bash

set -Eeu -o pipefail

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

TIMESTAMP=$(date "+%Y-%m-%d_%H-%M-%S")

MASSES=(
    100
    200
    300
    400
)

CTAUS=(
    3
    30
    50
    100
    300
)

PATTERN_MASS=$(printf "|MStau-%s" "${MASSES[@]}")
PATTERN_MASS=${PATTERN_MASS:1}

PATTERN_CTAU=$(printf "|ctau-%smm_mLSP" "${CTAUS[@]}")
PATTERN_CTAU=${PATTERN_CTAU:1}

INDIR=results/limits${SUFFIX}/llstau_${TYPE}/channels_BRT2/eras_all

for d in $(find "${INDIR}/" -mindepth 1 -maxdepth 1 -type d | grep -E "${PATTERN_MASS}" | grep -E "${PATTERN_CTAU}" | sort -V); do
    echo $d
    ./fit_limits.py --toylimits --indir $d $ARGS --timestamp $TIMESTAMP
done
