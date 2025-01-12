#!/usr/bin/env python3

import os

import ROOT
ROOT.gROOT.SetBatch(True)

import utils.utils


def add_hists(h1, h2, c1, c2) :
    
    htmp = h1.Clone()
    htmp.Add(h1, h2, c1, c2)
    
    #print(type(htmp))
    return htmp


def main() :
    
    infilename = "output/stau_kinematics_study/output_stau_kinematics_study.root"
    infile = ROOT.TFile.Open(infilename, "READ")
    
    outdir = os.path.dirname(infilename)
    os.system(f"mkdir -p {outdir}")
    
    l_stautype = [
        "stau_LH",
        "stau_RH",
        "stau_MM",
    ]
    
    l_stautype_comb = [
        #"stau_MD: stau_MD = {stau_LH}.Clone(); stau_MD.Add({stau_LH}, {stau_RH}, 0.00933, 0.003591)",
        "stau_MD: add_hists({stau_LH}, {stau_RH}, 0.00933, 0.003591)",
    ]
    
    l_histname = [
        #"GenVisTauh_pt",
        #"GenVisTaul_pt",
        #"GenVisTaul_e_by_stau_mass",
        
        "GenVisTauh_pt",
        "GenVisTauh_ebym",
        
        "GenVisTauh1_pt",
        "GenVisTauh2_pt",
        
        "Jet1_pt",
        "Jet2_pt",
        
        "MET_pt",
    ]
    
    d_xrange = {}
    #d_xrange["GenVisTauh_pt"] = (0, 500)
    #d_xrange["GenVisTaul_pt"] = (0, 500)
    #d_xrange["GenVisTaul_e_by_stau_mass"] = (0, 0.5)
    
    d_xrange["GenVisTauh_pt"] = (0, 500)
    d_xrange["GenVisTauh_ebym"] = (0, 1)
    d_xrange["GenVisTauh1_pt"] = (0, 500)
    d_xrange["GenVisTauh2_pt"] = (0, 500)
    d_xrange["Jet1_pt"] = (0, 500)
    d_xrange["Jet2_pt"] = (0, 500)
    d_xrange["MET_pt"] = (0, 500)
    
    
    d_yrange = {}
    #d_yrange["GenVisTauh_pt"] = (0, 0.1)
    #d_yrange["GenVisTaul_pt"] = (0, 0.2)
    #d_yrange["GenVisTaul_e_by_stau_mass"] = (0, 0.05)
    
    d_yrange["GenVisTauh_pt"] = (0, 0.1)
    d_yrange["GenVisTauh_ebym"] = (0, 0.05)
    d_yrange["GenVisTauh1_pt"] = (0, 0.1)
    d_yrange["GenVisTauh2_pt"] = (0, 0.2)
    d_yrange["Jet1_pt"] = (0, 0.1)
    d_yrange["Jet2_pt"] = (0, 0.2)
    d_yrange["MET_pt"] = (0, 0.1)
    
    
    d_hist = {}
    
    for stautype in l_stautype :
        
        d_hist[stautype] = {}
        
        for histname in l_histname :
            
            hist = infile.Get(f"{stautype}/{histname}").Clone()
            hist.SetDirectory(0)
            
            hist.Scale(1.0 / hist.Integral(0, hist.GetNbinsX()+1))
            print(hist.GetEntries(), hist.Integral())
            
            d_hist[stautype][histname] = hist
    
    for stautype in l_stautype_comb :
        
        stautype, operation = stautype.split(":")
        stautype = stautype.strip()
        operation = operation.strip()
        
        d_hist[stautype] = {}
        
        for histname in l_histname :
            
            operation_tmp = operation
            
            for key in l_stautype:
                
                operation_tmp = operation_tmp.replace(f"{{{key}}}", f"d_hist['{key}']['{histname}']")
            
            print(operation_tmp)
            hist = eval(operation_tmp)
            hist.Scale(1.0 / hist.Integral(0, hist.GetNbinsX()+1))
            hist.SetDirectory(0)
            
            print(type(hist))
            d_hist[stautype][histname] = hist
    
    infile.Close()
    
    for histname in l_histname :
        
        l_hist = []
        outfilename = f"{outdir}/{histname}.pdf"
        
        xrange = d_xrange[histname]
        yrange = d_yrange[histname]
        
        #for istautype, stautype in enumerate(l_stautype) :
        for istautype, stautype in enumerate(d_hist.keys()) :
            
            print(stautype, histname)
            hist = d_hist[stautype][histname].Clone()
            hist.SetTitle(stautype)
            hist.SetMarkerColor(istautype+1)
            hist.SetLineColor(istautype+1)
            hist.SetLineWidth(3)
            hist.SetMarkerSize(0)
            #hist.SetDrawOption("hist E2")
            hist.SetOption("hist")
            
            l_hist.append(hist)
        
        utils.utils.root_plot1D_legacy(
            l_hist = l_hist,
            #ratio_num_den_pairs = [(h1_num, h1_den)],
            outfile = outfilename,
            xrange = xrange,
            yrange = yrange,
            logx = False, logy = False,
            xtitle = histname,
            ytitle = "a.u.",
            #xtitle_ratio = d_hist[histName]["xtitle"], ytitle_ratio = "Efficiency",
            centertitlex = True, centertitley = True,
            centerlabelx = False, centerlabely = False,
            gridx = True, gridy = False,
            ndivisionsx = None,
            stackdrawopt = "nostack",
            legendpos = "UR",
            legendtitle = "[m_{#tilde{#tau}}(250), m_{LSP}(1)]",
            legendncol = 1,
            #legendtextsize = 0.04,
            legendwidthscale = 0.8,
            #legendheightscale = 1.5,
            lumiText = "2018 (13 TeV)"
        )


if __name__ == "__main__" :
    
    main()