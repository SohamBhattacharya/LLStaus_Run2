#!/usr/bin/env python3

import array
import numpy
import os

import ROOT
ROOT.gROOT.SetBatch(True)

import utils.utils

#era = "2016_preVFP"
era = "2016_postVFP"
#era = "2017"
#era = "2018"

fname_num = f"output/pepper_condor/MET-trigger-efficiency/{era}/hists/Cut-008_passed_probe_triggers_MET.root"
fname_den = f"output/pepper_condor/MET-trigger-efficiency/{era}/hists/Cut-007_has_2_jets_MET.root"

outdir = f"output/trigger_sf_study/{era}"
os.system(f"mkdir -p {outdir}")

infile_num = ROOT.TFile.Open(fname_num)
infile_den = ROOT.TFile.Open(fname_den)

l_histname_mc = [
    "WNJetsToLNu",
]

l_histname_data = [
    f"SingleMuon_{era}",
]

h1_num_mc = None
h1_num_data = None
h1_den_mc = None
h1_den_data = None

a_bins = array.array(
    "d",
    list(range(0, 100, 20)) + list(range(100, 300, 10)) + list(range(300, 400, 20)) + list(range(400, 600, 50)) + list(range(600, 1000, 200)) + [1000, 2000]
)

print(a_bins)

for hname in l_histname_mc :
    
    if(h1_num_mc is None) :
        h1_num_mc = infile_num.Get(hname).Clone("h1_num_mc")
        h1_num_mc.SetDirectory(0)
    else :
        h1_num_mc.Add(infile_num.Get(hname).Clone())
    
    if(h1_den_mc is None) :
        h1_den_mc = infile_den.Get(hname).Clone("h1_den_mc")
        h1_den_mc.SetDirectory(0)
    else :
        h1_den_mc.Add(infile_den.Get(hname).Clone())


for hname in l_histname_data :
    
    if(h1_num_data is None) :
        h1_num_data = infile_num.Get(hname).Clone("h1_num_data")
        h1_num_data.SetDirectory(0)
    else :
        h1_num_data.Add(infile_num.Get(hname).Clone())
    
    if(h1_den_data is None) :
        h1_den_data = infile_den.Get(hname).Clone("h1_den_data")
        h1_den_data.SetDirectory(0)
    else :
        h1_den_data.Add(infile_den.Get(hname).Clone())

d_hist = {
    "h1_num_mc": h1_num_mc,
    "h1_num_data": h1_num_data,
    "h1_den_mc": h1_den_mc,
    "h1_den_data": h1_den_data,
}

#h1_num_mc.Rebin(1, f"{h1_num_mc.GetName()}_rebinned", a_bins)
#h1_num_data.Rebin(1, f"{h1_num_data.GetName()}_rebinned", a_bins)
#h1_den_mc.Rebin(1, f"{h1_den_mc.GetName()}_rebinned", a_bins)
#h1_den_data.Rebin(1, f"{h1_den_data.GetName()}_rebinned", a_bins)

#h1_num_mc = h1_num_mc.Rebin(len(a_bins)-1, "h1_num_mc", a_bins)
#h1_num_data = h1_num_data.Rebin(len(a_bins)-1, "h1_num_data", a_bins)
#h1_den_mc = h1_den_mc.Rebin(len(a_bins)-1, "h1_den_mc", a_bins)
#h1_den_data = h1_den_data.Rebin(len(a_bins)-1, "h1_den_data", a_bins)

for key, hist in d_hist.items() :
    
    hist = hist.Rebin(len(a_bins)-1, "", a_bins)
    nBin = hist.GetNbinsX()
    
    # Set the last and overflow bins to the same value (sum of these two bins)
    val1 = hist.GetBinContent(nBin)
    val2 = hist.GetBinContent(nBin+1)
    hist.AddBinContent(nBin, val2)
    hist.AddBinContent(nBin+1, val1)
    
    d_hist[key] = hist

teff_mc = ROOT.TEfficiency(d_hist["h1_num_mc"], d_hist["h1_den_mc"])
teff_mc.SetName("teff_mc")
teff_data = ROOT.TEfficiency(d_hist["h1_num_data"], d_hist["h1_den_data"])
teff_data.SetName("teff_data")

h1_eff_mc = d_hist["h1_num_mc"].Clone("h1_eff_mc")
h1_eff_mc.Divide(d_hist["h1_num_mc"], d_hist["h1_den_mc"])#, 1, 1, "B")
h1_eff_mc.SetTitle("MC")
h1_eff_mc.SetMarkerColor(2)
h1_eff_mc.SetMarkerSize(2)
h1_eff_mc.SetMarkerStyle(20)
h1_eff_mc.SetLineColor(2)
h1_eff_mc.SetLineWidth(2)
h1_eff_mc.SetDrawOption("PE1")

h1_eff_data = d_hist["h1_num_data"].Clone("h1_eff_data")
h1_eff_data.Divide(d_hist["h1_num_data"], d_hist["h1_den_data"])#, 1, 1, "B")
h1_eff_data.SetTitle("Data")
h1_eff_data.SetMarkerColor(1)
h1_eff_data.SetMarkerSize(2)
h1_eff_data.SetMarkerStyle(20)
h1_eff_data.SetLineColor(1)
h1_eff_data.SetLineWidth(2)
h1_eff_data.SetDrawOption("PE1")

# Set the error on the efficiency to the binomial error from the TEfficiency object
nBin = h1_eff_mc.GetNbinsX()
for iBin in range(1, nBin+2) :
    
    #print(iBin, h1_sf.GetXaxis().GetBinLowEdge(iBin), h1_sf.GetXaxis().GetBinUpEdge(iBin), h1_sf.IsBinOverflow(iBin), h1_sf.GetBinContent(iBin))
    
    h1_eff_mc.SetBinError(iBin, 0.5*(teff_mc.GetEfficiencyErrorUp(iBin) + teff_mc.GetEfficiencyErrorLow(iBin)))
    h1_eff_data.SetBinError(iBin, 0.5*(teff_data.GetEfficiencyErrorUp(iBin) + teff_data.GetEfficiencyErrorLow(iBin)))

h1_sf = h1_eff_data.Clone("h1_sf")
h1_sf.Divide(h1_eff_mc)

fout = ROOT.TFile.Open(f"{outdir}/trigger-eff-sfs_MET.root", "RECREATE")
fout.cd()

# Write
for key, hist in d_hist.items() :
    
    hist.Write()

h1_eff_mc.Write()
h1_eff_data.Write()

h1_sf.Write()

teff_mc.Write()
teff_data.Write()

fout.Close()

utils.utils.root_plot1D_legacy(
    l_hist = [h1_eff_data, h1_eff_mc],
    ratio_num_den_pairs = [(h1_eff_data, h1_eff_mc)],
    outfile = f"{outdir}/trigger-eff-sfs_MET.pdf",
    xrange = (0, 800),
    yrange = (0, 1.2),
    logx = False, logy = False,
    ytitle = "Efficiency",
    xtitle_ratio = "p^{miss}_{T} [GeV]",
    ytitle_ratio = "Data / MC",
    yrange_ratio = (0.5, 1.5),
    centertitlex = True, centertitley = True,
    centerlabelx = False, centerlabely = False,
    gridx = True, gridy = True,
    ndivisionsx = None,
    ndivisionsy_ratio = (4, 5, 0), 
    stackdrawopt = "nostack",
    ratiodrawopt = "PE1",
    legendpos = "UL",
    legendncol = 1,
    legendwidthscale = 1.3,
    legendheightscale = 1.5,
    lumiText = f"{era} (13 TeV)",
)
