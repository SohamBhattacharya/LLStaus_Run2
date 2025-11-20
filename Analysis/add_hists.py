#!/usr/bin/env python3

import ROOT


def main() :
    
    #infname = "/home/soham/mnt/desy_dust/sobhatta/work/LongLivedStaus/shedprog_2017-studies/LLStaus_Run2/Analysis/output/2017/output_wjets/output_wjets_v2/hists/Cut_016_has_more_two_jets_jet_dxy.root"
    #l_hists = [
    #    "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/hist",
    #    "W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/hist",
    #    "W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/hist",
    #    "W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/hist",
    #    "W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/hist",
    #]
    #outhistname = "WNJetsToLNu"
    #outfname = "tmp/"+"output/2017/output_wjets/output_wjets_v2/hists/Cut_016_has_more_two_jets_jet_dxy.root".replace("/", "_")
    
    
    infname = "/home/soham/mnt/desy_dust/sobhatta/work/LongLivedStaus/shedprog_2017-studies/LLStaus_Run2/Analysis/output/2018/output_wjets/output_wjets_v2/hists/Cut_016_has_more_two_jets_jet_dxy.root"
    l_hists = [
        "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/hist",
        "W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/hist",
        "W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/hist",
        "W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/hist",
        "W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/hist",
    ]
    outhistname = "WNJetsToLNu"
    outfname = "tmp/"+"output/2018/output_wjets/output_wjets_v2/hists/Cut_016_has_more_two_jets_jet_dxy.root".replace("/", "_")
    
    
    infile = ROOT.TFile(infname, "READ")
    outfile = ROOT.TFile(outfname, "RECREATE")
    
    outhist = None
    
    for hname in l_hists :
        
        inhist = infile.Get(hname)
        inhist.SetDirectory(0)
        
        if outhist is None :
            outhist = inhist.Clone()
            outhist.SetName(outhistname)
            outhist.SetTitle(outhistname)
            outhist.SetDirectory(0)
        else :
            outhist.Add(inhist)
        
        outfile.WriteTObject(inhist, hname.replace("/", "_"))
    
    outhist_norm = outhist.Clone(f"{outhist.GetName()}_norm")
    outhist_norm.SetTitle(f"{outhist_norm.GetTitle()}")
    outhist_norm.Scale(1.0/outhist.GetEntries())
    outhist_norm.SetDirectory(0)
    
    outfile.cd()
    outhist.Write()
    outhist_norm.Write()
    outfile.Close()
    infile.Close()
    print(f"Output hist saved to {outfname}")
    
    return 0


if __name__ == "__main__":
    
    main()