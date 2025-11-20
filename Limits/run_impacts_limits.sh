
#!/bin/bash

# Exit if a command fails
set -Eeu -o pipefail

# Prints command before executing
set -o xtrace

DIR="$1"
#SIGNAL="$2" # 0 or 1

if [ $# -ge 2 ] && [ -n "$2" ]; then
    SIGNAL="$2"
else
    SIGNAL=""
fi

if [ -z "${DIR}" ]; then
    echo "Error. Usage:"
    echo "[ARG] is required, <ARG> is optional"
    echo "run_impacts_limits.sh [DIRECTORY] <SIGNAL: 0/1>"
fi

NPARALLEL=5

# --there does not work with combineTool.py -M Impacts
# Need to cd to the directory and run impacts there

WSPACE="workspace.root"

task(){
    indir=$1
    wd="${indir}/impacts"
    #args=""
    args="--rMin 0"
    #args="--rMin -10"
    
    if [ -n "${SIGNAL}" ]; then
        wd="${wd}_expectSignal${SIGNAL}"
        args="--expectSignal ${SIGNAL} -t -1 --rMin -100"
    fi
    
    mkdir -p $wd
    pushd $wd
    cp ../$WSPACE ./
    
    # First Stage: obtain the best fit for each POI with all nuisance profiling
    
    nice -n 10 combineTool.py \
    -M Impacts \
    -d $WSPACE \
    -m 120 \
    --doInitialFit \
    --robustFit 1 \
    --setRobustFitAlgo Minuit2 \
    --maxFailedSteps 100 \
    --forceRecreateNLL \
    --parallel 15 \
    $args

    # Second Stage: fit scan for each nuisance parameter

    nice -n 10 combineTool.py \
    -M Impacts \
    -d $WSPACE \
    -m 120 \
    --doFits \
    --robustFit 1 \
    --setRobustFitAlgo Minuit2 \
    --maxFailedSteps 100 \
    --forceRecreateNLL \
    --parallel 15 \
    $args

    # Collect outputs

    nice -n 10 combineTool.py \
    -M Impacts \
    -d $WSPACE \
    -m 120 \
    -o impacts.json \
    --parallel 15

    # Plot pulls and impacts

    nice -n 10 plotImpacts.py \
    -i impacts.json \
    -o impacts \
    --summary \
    --cms-label "Private work"
    
    popd
}

for wd in $(find $DIR -mindepth 0 -maxdepth 0 -type d); do
    
    task "$wd" &
    #if [[ $(jobs -r -p | wc -l) -ge $NPARALLEL ]]; then wait; fi
    
    while [[ $(jobs -r -p | wc -l) -ge $NPARALLEL ]]; do
        sleep 1
    done
    
done

wait

#--job-mode condor \
#--task-name $TASKNAME \
#--merge 1 \
#--sub-opts='+RequestRuntime = 172800\nRequestMemory = 31000' \
