#!/usr/bin/env -S python -u

import array
import numpy
import os
import re
import sortedcontainers
import ROOT

import utils.utils as utils
import utils.commonutils as cmut


def get_binning(start, end, step, include_end = True) :
    
    bins = numpy.arange(start, end, step)
    
    if (include_end) :
        
        if (bins[-1] != end) :
            
            bins = numpy.append(bins, [end])
        
        else :
            
            bins = bins[0: -1]
    
    else :
        
        if (bins[-1] == end) :
            
            bins = bins[0: -1]
    
    return bins


def main() :
    
    fname_in = "output/stau_samples_compare/output_stau_samples_compare.root"
    
    d_hist_config = sortedcontainers.SortedDict()
    
    d_hist_config["GenStau1_pt"] = {
        "label": "p_{T}(#tilde{#tau}_{1}) [GeV]",
        "rebin": get_binning(0, 1000, 50),
        "xrange": (0, 1000),
    }
    
    d_hist_config["GenStau2_pt"] = {
        "label": "p_{T}(#tilde{#tau}_{2}) [GeV]",
        "rebin": get_binning(0, 1000, 50),
        "xrange": (0, 1000),
    }
    
    d_hist_config["GenStau1_pathL"] = {
        "label": "l(#tilde{#tau}_{1}) [cm]",
        "rebin": get_binning(0, 10, 0.1),
        "xrange": (0, 10),
    }
    
    d_hist_config["GenStau2_pathL"] = {
        "label": "l(#tilde{#tau}_{2}) [cm]",
        "rebin": get_binning(0, 10, 0.1),
        "xrange": (0, 10),
    }
    
    d_hist_config["GenTau1_pt"] = {
        "label": "p_{T}(#tau_{1}) [GeV]",
        "rebin": get_binning(0, 1000, 50),
        "xrange": (0, 1000),
    }
    
    d_hist_config["GenTau2_pt"] = {
        "label": "p_{T}(#tau_{2}) [GeV]",
        "rebin": get_binning(0, 1000, 50),
        "xrange": (0, 1000),
    }
    
    d_hist_config["GenTau1_vertexR"] = {
        "label": "R_{vtx}(#tau_{1}) [cm]",
        "rebin": get_binning(0, 10, 0.1),
        "xrange": (0, 10),
    }
    
    d_hist_config["GenTau2_vertexR"] = {
        "label": "R_{vtx}(#tau_{2}) [cm]",
        "rebin": get_binning(0, 10, 0.1),
        "xrange": (0, 10),
    }
    
    d_hist_config["Jet1_pt"] = {
        "label": "p_{T}(j1) [GeV]",
        "rebin": get_binning(0, 1000, 50),
        "xrange": (0, 1000),
    }
    
    d_hist_config["Jet2_pt"] = {
        "label": "p_{T}(j2) [GeV]",
        "rebin": get_binning(0, 1000, 50),
        "xrange": (0, 1000),
    }
    
    d_hist_config["Jet1_dxy"] = {
        "label": "d_{xy}(j1) [cm]",
        "rebin": get_binning(0, 5, 0.1),
        "xrange": (0, 5),
    }
    
    d_hist_config["Jet2_dxy"] = {
        "label": "d_{xy}(j2) [cm]",
        "rebin": get_binning(0, 5, 0.1),
        "xrange": (0, 5),
    }
    
    d_hist_config["Jet1_disTauTag_score1"] = {
        "label": "DisTau score (j1)",
        "rebin": get_binning(0, 1, 0.1),
        "xrange": (0, 1),
    }
    
    d_hist_config["Jet2_disTauTag_score1"] = {
        "label": "DisTau score (j2)",
        "rebin": get_binning(0, 1, 0.1),
        "xrange": (0, 1),
    }
    
    d_hist_config["MET_pt"] = {
        "label": "p^{miss}_{T} [GeV]",
        "rebin": get_binning(0, 1000, 50),
        "xrange": (0, 1000),
    }
    
    
    d_stau_config = sortedcontainers.SortedDict()
    
    d_stau_config["MStau-250_ctau-100mm"] = [
        "stau250_lsp1_ctau100mm_ul18_private",
        "stau250_lsp1_ctau100mm_ul18_central",
    ]
    
    infile = ROOT.TFile.Open(fname_in)
    outdir = os.path.dirname(fname_in)
    
    if (outdir) :
        
        os.system(f"mkdir -p {outdir}")
    
    for hist_key, hist_config in d_hist_config.items() :
        
        for stau_key, stau_config in d_stau_config.items() :
            
            l_hist = []
            
            for cfg_name in stau_config :
                
                hist_name = f"{cfg_name}/{hist_key}"
                print(hist_name)
                
                hist_tmp = cmut.get_hist(
                    histfile = infile,
                    histname = hist_key,
                    samples = [cfg_name],
                    rebin = hist_config["rebin"],
                )
                
                hist_tmp.Scale(1.0 / hist_tmp.Integral())
                
                #hist_tmp = infile.Get(hist_name)
                #hist_tmp.SetDirectory(0)
                
                color = utils.ColorIterator(len(l_hist), 0)
                
                #regexp = "MStau-(?P<mstau>\d+)_ctau-(?P<ctau>\w+)mm"
                regexp = "stau(?P<mstau>\d+)_lsp1_ctau(?P<ctau>\w+)mm_ul18_(?P<sample>\w+)"
                rgx = re.compile(regexp)
                
                stau_info = cmut.parse_string_regex(
                    s = cfg_name,
                    regexp = regexp,
                )
                #stau_info = [m.groupdict() for m in rgx.finditer(cfg_name)][0]
                
                mstau = stau_info["mstau"]
                ctau = float(stau_info["ctau"].replace("p", "."))
                sample = stau_info["sample"]
                
                hist_title = f"#tilde{{#tau}}_{{1}}({mstau}), c#tau_{{0}}({ctau:0.1f} mm) [{sample}]"
                
                if "from" in cfg_name :
                    
                    ctaufrom = 0.1*float(stau_info["ctaufrom"].replace("p", "."))
                    hist_title = f"{hist_title}#Leftarrowc#tau_{{0}}({ctaufrom:0.1f} mm)"
                
                hist_tmp.SetTitle(hist_title)
                hist_tmp.SetLineWidth(2)
                hist_tmp.SetLineColor(color)
                hist_tmp.SetMarkerSize(0)
                hist_tmp.SetMarkerColor(color)
                hist_tmp.SetFillStyle(0)
                hist_tmp.SetOption("hist E")
                l_hist.append(hist_tmp)
            
            l_ratio_pairs = [(l_hist[_i], l_hist[0]) for _i in range(1, len(l_hist))]
            
            outfile = f"{outdir}/{stau_key}/{hist_key}.pdf"
            
            utils.root_plot1D_legacy(
                l_hist = l_hist,
                ratio_num_den_pairs = l_ratio_pairs,
                outfile = outfile,
                xrange = hist_config["xrange"],
                yrange = (1e-4, 10),
                logx = False,
                logy = True,
                ytitle = "a.u.",
                xtitle_ratio = hist_config["label"],
                ytitle_ratio = "Ratio",
                yrange_ratio = (0, 2),
                ratiodrawopt = "hist E",
                centertitlex = True, centertitley = True,
                centerlabelx = False, centerlabely = False,
                gridx = True, gridy = True,
                ndivisionsx = None,
                ndivisionsy_ratio = (4, 5, 0),
                stackdrawopt = "nostack",
                legendpos = "UR",
                #legendtitle = "#splitline{%s}{(%s)}" %(dataset, tag),
                legendncol = 1,
                legendtextsize = 0.04,
                legendwidthscale = 1.2,
                legendheightscale = 1.5,
                lumiText = "2018 (13 TeV)",
            )
    
    infile.Close()
    
    return 0

if __name__ == "__main__" :
    
    main()