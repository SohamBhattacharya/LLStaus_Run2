#!/bin/bash

NAME=$1

#echo "preVFP:"
#dasgoclient -query="dataset=/$1/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11*/MINIAODSIM" | sort -V
#
#echo ""
#
#echo "postVFP:"
#dasgoclient -query="dataset=/$1/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17*/MINIAODSIM" | sort -V


dasgoclient -query="dataset=/$1/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9*/MINIAODSIM" | sort -V