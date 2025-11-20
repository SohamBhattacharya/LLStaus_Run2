#!/usr/bin/env python3

import os
import json

import ROOT
ROOT.gROOT.SetBatch(True)

import utils.utils as utils
import utils.commonutils as cmut


def main() :
    
    infilename = "output/stau_score_study/output_stau_score_study.root"
    #infilename = "output/stau_score_study_dxy-gt-0p2/output_stau_score_study.root"
    infile = ROOT.TFile.Open(infilename, "READ")
    
    outdir = os.path.dirname(infilename)
    os.system(f"mkdir -p {outdir}")
    
    d_sample_xsecs = {}
    
    with open("configs/crosssections.json", "r") as fopen:
        d_sample_xsecs = json.load(fopen)
    
    d_processes = {}
    d_processes["QCD"] = [
        "QCD_Pt_30to50_TuneCP5_13TeV_pythia8",
        "QCD_Pt_50to80_TuneCP5_13TeV_pythia8",
        "QCD_Pt_80to120_TuneCP5_13TeV_pythia8",
        "QCD_Pt_120to170_TuneCP5_13TeV_pythia8",
        "QCD_Pt_170to300_TuneCP5_13TeV_pythia8",
        "QCD_Pt_300to470_TuneCP5_13TeV_pythia8",
        "QCD_Pt_470to600_TuneCP5_13TeV_pythia8",
        "QCD_Pt_600to800_TuneCP5_13TeV_pythia8",
        "QCD_Pt_800to1000_TuneCP5_13TeV_pythia8",
        "QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8",
        "QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8",
        "QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8",
        "QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8",
        "QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8",
    ]
    
    l_plot_processes = [
        "TTToHadronic",
        #"TTToHadronic:light",
        #"TTToHadronic:b",
        #"TTToHadronic:c",
        #"TTToHadronic:bc",
        
        #"QCD:light",
        #"QCD:bc",
        
        "MStau-200_ctau-5mm",
        #"MStau-200_ctau-50mm",
        "MStau-200_ctau-100mm",
        #"MStau-200_ctau-1000mm",
        
        "MStau-400_ctau-5mm",
        "MStau-400_ctau-100mm",
        #"MStau-400_ctau-1000mm",
    ]
    
    d_proc_label = {}
    d_proc_label["TTToHadronic"] = "t#bar{t}"
    d_proc_label["TTToHadronic:light"] = "t#bar{t} (u/d/s/g jets)"
    d_proc_label["TTToHadronic:b"] = "t#bar{t} (b jets)"
    d_proc_label["TTToHadronic:c"] = "t#bar{t} (c jets)"
    d_proc_label["TTToHadronic:bc"] = "t#bar{t} (b/c jets)"
    
    d_proc_label["QCD:light"] = "Multijet (u/d/s/g jets)"
    d_proc_label["QCD:bc"] = "Multijet (b/c jets)"
    
    d_proc_label["MStau-200_ctau-5mm"] = "m_{#tilde{#tau}}=200 GeV, c#tau_{0}=5 mm"
    d_proc_label["MStau-200_ctau-50mm"] = "m_{#tilde{#tau}}=200 GeV, c#tau_{0}=50 mm"
    d_proc_label["MStau-200_ctau-100mm"] = "m_{#tilde{#tau}}=200 GeV, c#tau_{0}=100 mm"
    d_proc_label["MStau-200_ctau-1000mm"] = "m_{#tilde{#tau}}=200 GeV, c#tau_{0}=1000 mm"
    
    d_proc_label["MStau-400_ctau-5mm"] = "m_{#tilde{#tau}}=400 GeV, c#tau_{0}=5 mm"
    d_proc_label["MStau-400_ctau-100mm"] = "m_{#tilde{#tau}}=400 GeV, c#tau_{0}=100 mm"
    d_proc_label["MStau-400_ctau-1000mm"] = "m_{#tilde{#tau}}=400 GeV, c#tau_{0}=1000 mm"
    
    
    d_proc_color = {}
    d_proc_color["TTToHadronic"] = utils.get_cms_colors(0)
    d_proc_color["TTToHadronic:light"] = utils.get_cms_colors(0)
    d_proc_color["TTToHadronic:b"] = utils.get_cms_colors(0)
    d_proc_color["TTToHadronic:c"] = utils.get_cms_colors(0)
    d_proc_color["TTToHadronic:bc"] = utils.get_cms_colors(0)
    
    d_proc_color["QCD:light"] = utils.get_cms_colors(1)
    d_proc_color["QCD:bc"] = utils.get_cms_colors(1)
    
    d_proc_color["MStau-200_ctau-5mm"] = utils.get_cms_colors(2)
    d_proc_color["MStau-200_ctau-100mm"] = utils.get_cms_colors(2)
    
    d_proc_color["MStau-400_ctau-5mm"] = utils.get_cms_colors(4)
    d_proc_color["MStau-400_ctau-100mm"] = utils.get_cms_colors(4)
    
    
    d_proc_linestyle = {}
    d_proc_linestyle["TTToHadronic"] = 1
    d_proc_linestyle["TTToHadronic:light"] = 1
    d_proc_linestyle["TTToHadronic:b"] = 2
    d_proc_linestyle["TTToHadronic:c"] = 3
    d_proc_linestyle["TTToHadronic:bc"] = 7
    
    d_proc_linestyle["QCD:light"] = 1
    d_proc_linestyle["QCD:bc"] = 7

    d_proc_linestyle["MStau-200_ctau-5mm"] = 1
    d_proc_linestyle["MStau-200_ctau-100mm"] = 7

    d_proc_linestyle["MStau-400_ctau-5mm"] = 1
    d_proc_linestyle["MStau-400_ctau-100mm"] = 7
    
    histname_nevents = "nEvents"
    
    d_histname = {}
    d_histname["Jet_disTauTag_score1"] = {
        "TTToHadronic": "Jet_light_disTauTag_score1",
        "TTToHadronic:light": "Jet_light_disTauTag_score1",
        "TTToHadronic:b": "Jet_b_disTauTag_score1",
        "TTToHadronic:c": "Jet_c_disTauTag_score1",
        "TTToHadronic:bc": "Jet_bc_disTauTag_score1",
        
        "QCD:light": "Jet_light_disTauTag_score1",
        "QCD:bc": "Jet_bc_disTauTag_score1",
        
        "MStau-200_ctau-5mm": "Jet_disTauTag_score1",
        "MStau-200_ctau-100mm": "Jet_disTauTag_score1",
        
        "MStau-400_ctau-5mm": "Jet_disTauTag_score1",
        "MStau-400_ctau-100mm": "Jet_disTauTag_score1",
    }
    
    d_xrange = {}
    d_xrange["Jet_disTauTag_score1"] = (0, 1)
    
    d_yrange = {}
    d_yrange["Jet_disTauTag_score1"] = (1e-5, 1e2)
    
    d_xtitle = {}
    d_xtitle["Jet_disTauTag_score1"] = "Displaced#scale[0.4]{ }#tau_{h} tagger score"
    
    d_hist = {}
    
    rebin = cmut.get_binning(0, 1, 0.05)
    
    for proc in l_plot_processes :
        
        d_hist[proc] = {}
        procname = proc.split(":")[0]
        
        l_samples = d_processes.get(procname, [procname])
        
        for histkey in list(d_histname.keys()) :
            
            for sample in l_samples :
                
                histname = d_histname[histkey][proc]
                print(f"[{sample}] [{histkey}] [{histname}]")
                #print(f"[{sample}] [{histname}]")
                
                hist = cmut.get_hist(
                    histfile = infile,
                    histname = histname,
                    samples = [sample],
                    scales = [],
                    rebin = rebin,
                    underflow = True,
                    overflow = True,
                    set_min = None,
                    set_max = None,
                )
                
                hist_nevents = cmut.get_hist(
                    histfile = infile,
                    histname = histname_nevents,
                    samples = [sample],
                    scales = [],
                    rebin = None,
                    underflow = True,
                    overflow = True,
                    set_min = None,
                    set_max = None,
                )
                
                xsec = d_sample_xsecs.get(sample, 1)
                
                hist.Scale(xsec / hist_nevents.Integral())
                
                if histkey not in d_hist[proc] :
                    d_hist[proc][histkey] = hist.Clone()
                else :
                    d_hist[proc][histkey].Add(hist)
                
                bin_start = d_hist[proc][histkey].GetXaxis().FindBin(d_xrange[histkey][0])
                bin_end = d_hist[proc][histkey].GetXaxis().FindBin(d_xrange[histkey][1])
            
            d_hist[proc][histkey].Scale(1.0 / d_hist[proc][histkey].Integral(bin_start, bin_end))
        
        #for histkey in list(d_histname.keys()) :
        #    
        #    histname = d_histname[histkey][sample]
        #    print(f"[{sample}] [{histkey}] [{histname}]")
        #    #hist = infile.Get(f"{sample}/{histname}").Clone()
        #    #hist.SetDirectory(0)
        #    
        #    #hist = hist.Rebin(len(rebin)-1, "", rebin)
        #    
        #    hist = cmut.get_hist(
        #        histfile = infile,
        #        histname = histname,
        #        samples = samplename,
        #        scales = [],
        #        rebin = rebin,
        #        underflow = True,
        #        overflow = True,
        #        set_min = None,
        #        set_max = None,
        #    )
        #    
        #    print(hist.GetBinContent(0), hist.GetBinContent(hist.GetNbinsX()+1))
        #    #hist.Scale(1.0 / hist.Integral(0, hist.GetNbinsX()+1))
        #    print(hist.GetEntries(), hist.Integral())
        #    hist.Scale(1.0 / hist.Integral())
        #    print(hist.GetEntries(), hist.Integral())
        #    #print(hist.GetBinContent(5), hist.GetBinError(5), hist.GetBinError(5)/hist.GetBinContent(5))
        #    
        #    d_hist[sample][histkey] = hist
    
    for histkey in list(d_histname.keys()) :
        
        l_hist = []
        outfilename = f"{outdir}/{histkey}.pdf"
        
        xrange = d_xrange[histkey]
        yrange = d_yrange[histkey]
        
        for iproc, proc in enumerate(d_hist.keys()) :
            
            #color = utils.get_cms_colors(iproc)
            color = d_proc_color[proc]
            linestyle = d_proc_linestyle[proc]

            print(proc, histkey)
            hist = d_hist[proc][histkey].Clone()
            hist.SetTitle(d_proc_label[proc])
            hist.SetMarkerColor(color)
            hist.SetLineColor(color)
            hist.SetLineStyle(linestyle)
            hist.SetLineWidth(3)
            #hist.SetMarkerSize(0)
            hist.SetFillStyle(0)
            hist.SetDrawOption("hist")
            hist.SetOption("hist")
            
            l_hist.append(hist)
        
        utils.root_plot1D_legacy(
            l_hist = l_hist,
            #ratio_num_den_pairs = [(_hist, l_hist[-1]) for _hist in l_hist[:-1]],
            outfile = outfilename,
            xrange = xrange,
            yrange = yrange,
            logx = False,
            logy = True,
            xtitle = d_xtitle[histkey],
            ytitle = "a.u.",
            xtitle_ratio = d_xtitle[histkey],
            ytitle_ratio = f"Ratio wrt {l_hist[-1].GetTitle()}",
            centertitlex = True, centertitley = True,
            centerlabelx = False, centerlabely = False,
            #gridx = True, gridy = True,
            ndivisionsx = None,
            stackdrawopt = "nostack",
            ratiodrawopt = "E1",
            ndivisionsy_ratio = (2, 5, 0),
            legendpos = "UL",
            #legendtitle = "[m_{#tilde{#tau}}(250), m_{LSP}(1)]",
            legendncol = 1,
            legendtextsize = 0.045,
            legendwidthscale = 1.9,
            #legendheightscale = 0.5,
            legendheightscale = 1.2,
            legendpadleft_extra = 0.05,
            CMSextraText = "#kern[-0.15]{Simulation}",
            #lumiText = "2018 (13 TeV)"
            lumiText = "13 TeV"
        )


if __name__ == "__main__" :
    
    main()