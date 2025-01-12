#!/bin/bash

# 2016_preVFP
#TAG="RunIISummer20UL16MiniAODAPVv2-SUS_106X_mcRun2_asymptotic_preVFP_v11"

# 2016_postVFP
#TAG="RunIISummer20UL16MiniAODv2-SUS_106X_mcRun2_asymptotic_v17"

# 2017
TAG="RunIISummer20UL17MiniAODv2-SUS_106X_mc2017_realistic_v9"

# 2018
#TAG="RunIISummer20UL18MiniAODv2-SUS_106X_upgrade2018_realistic_v16_L1v1"

QUERY="dataset=/SMS-TStauStau_MStau*/$TAG*/MINIAODSIM"

for f in $(dasgoclient -query=`echo $QUERY` | grep -v Pilot | sort -V); do
    s=`echo $f | awk -F/ '{print $2}' | awk -F"_Tune" '{print $1}' | sed "s/-/_/g"`
    printf "%-80s%s\n" $(echo $s $f)
done
