#!/bin/bash

DIR="$1"

NPARALLEL=15

task(){
    wd=$1
    
    pushd $wd
    
    echo "Processing: $wd ..."
    
    nice -n 10 hadd -f -k higgsCombine.Test.HybridNew.mH120.merged.root `ls | grep -E "higgsCombine.Test.POINT.*.HybridNew.mH120.*[0-9]+.root" | sort -V`
    
    nice -n 10 combine workspace.root -M HybridNew --LHCmode LHC-limits --readHybridResults --grid=higgsCombine.Test.HybridNew.mH120.merged.root --plot=limit_scan_obs.png
    nice -n 10 combine workspace.root -M HybridNew --LHCmode LHC-limits --readHybridResults --grid=higgsCombine.Test.HybridNew.mH120.merged.root --plot=limit_scan_exp50.png --expectedFromGrid=0.5
    nice -n 10 combine workspace.root -M HybridNew --LHCmode LHC-limits --readHybridResults --grid=higgsCombine.Test.HybridNew.mH120.merged.root --plot=limit_scan_exp16.png --expectedFromGrid=0.16
    nice -n 10 combine workspace.root -M HybridNew --LHCmode LHC-limits --readHybridResults --grid=higgsCombine.Test.HybridNew.mH120.merged.root --plot=limit_scan_exp84.png --expectedFromGrid=0.84
    nice -n 10 combine workspace.root -M HybridNew --LHCmode LHC-limits --readHybridResults --grid=higgsCombine.Test.HybridNew.mH120.merged.root --plot=limit_scan_exp2p5.png --expectedFromGrid=0.025
    nice -n 10 combine workspace.root -M HybridNew --LHCmode LHC-limits --readHybridResults --grid=higgsCombine.Test.HybridNew.mH120.merged.root --plot=limit_scan_exp97p5.png --expectedFromGrid=0.975
    
    popd
}

for wd in $(find $DIR -mindepth 0 -maxdepth 0 -type d); do
    
    task "$wd" &
    
    while [[ $(jobs -r -p | wc -l) -ge $NPARALLEL ]]; do
        sleep 1
    done
    
done

wait