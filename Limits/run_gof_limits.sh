
#!/bin/bash

# Exit if a command fails
set -Eeu -o pipefail

# Prints command before executing
set -o xtrace

DIR="$1"
SIGNAL="$2"

#if [ $# -ge 2 ] && [ -n "$2" ]; then
#    SIGNAL="$2"
#else
#    SIGNAL=""
#fi

if [ -z "${DIR}" ] || [ -z "${SIGNAL}" ]; then
    echo "Error. Usage:"
    #echo "[ARG] is required, <ARG> is optional"
    echo "run_gof_limits.sh [DIRECTORY] [SIGNAL: b/sb/number]"
    #echo "Will run on data if no SIGNAL is given"
    echo "b = background-only, sb = signal+background, number = will set signal strength to this number"
fi

NPARALLEL=1

WSPACE="workspace.root"

task(){
    indir=$1
    wd="${indir}/gof"
    args=""
    args_toys=""
    title=""
    
    if [ "${SIGNAL}" = "b" ]; then
        wd="${wd}_${SIGNAL}"
        title="bkg-only"
        args="--freezeParameters r --setParameters r=0"
        args_toys="--freezeParameters r --setParameters r=0 --toysFrequentist -t 1 -s=-1:1000 --parallel 20"
    elif [ "${SIGNAL}" = "sb" ]; then
        wd="${wd}_${SIGNAL}"
        title="sig+bkg"
        args=""
        args_toys="--toysFrequentist -t 1 -s=-1:1000 --parallel 20"
    else
        suffix=$(sed "s/\./p/g" <<< "${SIGNAL}")
        wd="${wd}_s${suffix}"
        title="r=${SIGNAL}"
        args="--freezeParameters r --setParameters r=${SIGNAL}"
        args_toys="--freezeParameters r --setParameters r=${SIGNAL} --toysFrequentist -t 1 -s=-1:1000 --parallel 20"
    fi
    
    mkdir -p $wd
    pushd $wd
    cp ../$WSPACE ./
    
    # On data
    nice -n 10 combineTool.py \
    -M GoodnessOfFit \
    -d $WSPACE \
    --algo=saturated \
    --there \
    -n ".data" \
    $args
    
    # With toys
    nice -n 10 combineTool.py \
    -M GoodnessOfFit \
    -d $WSPACE \
    --algo=saturated \
    --there \
    -n ".toys" \
    $args_toys
    
    files=$(find higgsCombine.toys.GoodnessOfFit.mH120.*.root | sort -V)
    hadd -f higgsCombine.toys.GoodnessOfFit.root $files && rm $files
    
    nice -n 10 combineTool.py \
    -M CollectGoodnessOfFit \
    --input higgsCombine.data.GoodnessOfFit.mH120.root higgsCombine.toys.GoodnessOfFit.root \
    -o gof.json
    
    plotGof.py gof.json \
    --statistic saturated \
    --mass 120.0 \
    -o "gof_plot" \
    --title-right="${title}"
    
    popd
}

for wd in $(find $DIR -mindepth 0 -maxdepth 0 -type d); do
    
    task "$wd" &
    
    while [[ $(jobs -r -p | wc -l) -ge $NPARALLEL ]]; do
        sleep 1
    done
    
done

wait

#--job-mode condor \
#--task-name $TASKNAME \
#--merge 1 \
#--sub-opts='+RequestRuntime = 172800\nRequestMemory = 31000' \
