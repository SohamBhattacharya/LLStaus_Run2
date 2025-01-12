
#!/bin/bash

# Exit if a command fails
set -Eeu -o pipefail

# Prints command before executing
set -o xtrace

DIR="$1"
#DIR="tmp/test_DisTauSF_mass/DisTauSF/channels_all/eras_all/ZMT_wp-p*_dxy-gt-0p*"
#DIR="tmp/test_DisTauSF_mass/DisTauSF/channels_all/eras_all/ZMT_wp-p80_dxy-gt-0p00"
#DIR="tmp/test_DisTauSF_mass/DisTauSF/channels_all/eras_all/ZMT_wp-p95_dxy-gt-0p05"
#DIR="tmp/test_DisTauSF_mass/DisTauSF/channels_all/eras_all/ZMT_wp-p99_dxy-gt-0p07"

NPARALLEL=5

# --there does not work with combineTool.py -M Impacts
# Need to cd to the directory and run impacts there

WSPACE="workspace.root"

task(){
    wd=$1
    
    pushd $wd
    
    # First Stage: obtain the best fit for each POI with all nuisance profiling
    
    nice -n 10 combineTool.py \
    -M Impacts \
    -d $WSPACE \
    --mass 90 \
    --doInitialFit \
    --robustFit 1 \
    --cminDefaultMinimizerStrategy 0 \
    --redefineSignalPOIs SF \
    --setParameterRanges SF=0,3 \
    --forceRecreateNLL \
    -n ".DisTauSF" \
    --parallel 15

    # Second Stage: fit scan for each nuisance parameter

    nice -n 10 combineTool.py \
    -M Impacts \
    -d $WSPACE \
    --mass 90 \
    --doFits \
    --robustFit 1 \
    --cminDefaultMinimizerStrategy 0 \
    --redefineSignalPOIs SF \
    --setParameterRanges SF=0,3 \
    --forceRecreateNLL \
    -n ".DisTauSF" \
    --parallel 15

    # Collect outputs

    nice -n 10 combineTool.py \
    -M Impacts \
    -d $WSPACE \
    --mass 90 \
    -o impacts.json \
    -n ".DisTauSF" \
    --parallel 15

    # Plot pulls and impacts

    nice -n 10 plotImpacts.py \
    -i impacts.json \
    -o impacts \
    --POI SF \
    
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
