#!/usr/bin/env python3

import functools
import operator

import utils.commonutils as cmut

def main() :
    
    #cutflows_file = "/afs/desy.de/user/m/mykytaua/nfscms/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_6/2018/output_zmutau/zmutau_v1/cutflows.json"
    cutflows_file = "/home/soham/mnt/desy_dust/mykytaua/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_6/2018/output_zmutau/zmutau_v1/cutflows.json"
    
    xsec_file = "../Analysis/configs/crosssections.json"
    
    neventkey = "all.BeforeCuts"
    
    l_samples = [
        #"DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8",
        #"DY1JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8",
        #"DY2JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8",
        #"DY3JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8",
        #"DY4JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8",
        
        "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
        "W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
        "W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
        "W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
        "W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
    ]
    
    l_nevents_sample = []
    d_xsec = cmut.load_config(xsec_file)
    d_cutflows = cmut.load_config(cutflows_file)
    
    for sample in l_samples :
        
        #xsec = d_xsec[sample]
        nevents_sample = functools.reduce(operator.getitem, [sample]+neventkey.split("."), d_cutflows)
        
        cmut.logger.info(f"[{sample}] [{neventkey}] {nevents_sample} events")
        
        l_nevents_sample.append(nevents_sample)
    
    nevents_total = sum(l_nevents_sample)
    
    print(f"Total = {nevents_total}")
    
    return 0


if __name__ == "__main__" :
    
    main()