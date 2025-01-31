#!/bin/bash

DIR="$1"

NPARALLEL=5

task(){
    wd=$1
    
    pushd $wd
    
    echo "Processing: $wd ..."
    
    PostFitShapesFromWorkspace \
    -w workspace.root \
    --output postfit_s.root \
    --sampling \
    -f fitDiagnostics.limits.root:fit_s \
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