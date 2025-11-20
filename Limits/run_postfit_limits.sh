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
    -f fitDiagnostics.Test.root:fit_s \
    --postfit \
    --covariance
    
    PostFitShapesFromWorkspace \
    -w workspace.root \
    --output postfit_b.root \
    --sampling \
    -f fitDiagnostics.Test.root:fit_b \
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