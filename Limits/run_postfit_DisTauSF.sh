#!/bin/bash

DIR="$1"
#DIR="tmp/test_DisTauSF_mass_nbins1-60-80/DisTauSF/channels_all/eras_all/ZMT_wp-p*_dxy-gt-0p*"
#DIR="tmp/test_DisTauSF_mass_nbins2-60-70-80/DisTauSF/channels_all/eras_all/ZMT_wp-p*_dxy-gt-0p*"

NPARALLEL=5

task(){
    wd=$1
    
    pushd $wd
    
    echo "Processing: $wd ..."
    
    PostFitShapesFromWorkspace \
    -w workspace.root \
    --output postfit_s.root \
    --sampling \
    -f fitDiagnostics.DisTauSF.root:fit_s \
    --postfit \
    --covariance
    
    popd
}

for wd in $(find $DIR -mindepth 0 -maxdepth 0 -type d); do
    
    task "$wd" &
    
    while [[ $(jobs -r -p | wc -l) -ge $NPARALLEL ]]; do
        sleep 1
    done
    
done

wait