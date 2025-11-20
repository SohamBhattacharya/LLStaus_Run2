#!/usr/bin/env python3

import argparse
import numpy
import os
import sortedcontainers
import ROOT

import utils.utils as utils
import utils.commonutils as cmut

def main() :
    
    # Argument parser
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        "--dir",
        help = " ".join([
            "Input directory containing SMS-TStauStau_MStau-*.",
            "For e.g. results/limits_signal_v12/llstau_maximally-mixed/channels_BRT2/eras_all/",
        ]),
        type = str,
        required = True,
    )
    
    parser.add_argument(
        "--toystamp",
        help = "Timestamp (usually YYYY-MM-DD_hh-mm-ss) of the toy condor jobs",
        type = str,
        required = True,
    )
    
    parser.add_argument(
        "--limits",
        help = "Limits to compare",
        type = str,
        nargs = "+",
        required = False,
        default = [
            "obs",
            "exp",
            "exp_m1",
            "exp_p1",
            "exp_m2",
            "exp_p2",
        ]
    )
    
    
    # Parse arguments
    args = parser.parse_args()
    
    l_ctaus = numpy.loadtxt("configs/limits/sig_ctau_list.txt", dtype = str)
    #l_ctaus = [_val for _val in l_ctaus if _val in ["300mm"]]
    
    scenario = ""
    if "maximally-mixed" in args.dir :
        scenario = "Maximally mixed scenario"
    elif "mass-degenerate" in args.dir :
        scenario = "Mass degenerate scenario"
    
    infname_asymptotic =f"{args.dir}/limits/limits.root"
    infname_toys = f"{args.dir}/limits_toys/{args.toystamp}/limits.root"

    cmut.logger.info(f"Opening asymptotic limits file: {infname_asymptotic}")
    cmut.logger.info(f"Opening toy limits file: {infname_toys}")
    
    infile_asymptotic = ROOT.TFile.Open(infname_asymptotic, "READ")
    infile_toys = ROOT.TFile.Open(infname_toys, "READ")
    
    outdir = f"{args.dir}/compare_limits_asymptotic_toys"
    os.system(f"mkdir -p {outdir}")
    
    xmin, xmax = (90, 600)
    
    for ctau in l_ctaus :
        
        for limtype in args.limits :
            
            limtype_label = limtype
            limtype_label = limtype_label.replace("_m1", " #minus 1#sigma")
            limtype_label = limtype_label.replace("_p1", " #plus 1#sigma")
            limtype_label = limtype_label.replace("_m2", " #minus 2#sigma")
            limtype_label = limtype_label.replace("_p2", " #plus 2#sigma")

            gr_name = f"{ctau}/g1_xsecul_{limtype}_{ctau}"
            
            cmut.logger.info(f"Processing {gr_name} ...")
            
            g1_asymptotic = infile_asymptotic.Get(gr_name)
            g1_toys = infile_toys.Get(gr_name)
            
            g1_asymptotic_pm1 = None
            g1_toys_pm1 = None
            
            if limtype == "exp" :
                g1_asymptotic_pm1 = infile_asymptotic.Get(f"{ctau}/g1_xsecul_exp_pm1_{ctau}")
                g1_asymptotic_pm1.SetTitle(f"Asymptotic (exp #pm 1#sigma)")
                g1_asymptotic_pm1.SetFillColorAlpha(2, 0.3)
                #g1_asymptotic_pm1.SetFillStyle(3244)
                #g1_asymptotic_pm1.SetLineColor(2)
                #g1_asymptotic_pm1.SetLineWidth(2)
                g1_asymptotic_pm1.SetLineWidth(0)
                g1_asymptotic_pm1.SetOption("3")
                g1_asymptotic_pm1.GetHistogram().SetOption("3")
                g1_asymptotic_pm1.GetHistogram().SetOption("3")
                
                g1_toys_pm1 = infile_toys.Get(f"{ctau}/g1_xsecul_exp_pm1_{ctau}")
                g1_toys_pm1.SetTitle(f"Toys (exp #pm 1#sigma)")
                g1_toys_pm1.SetFillColorAlpha(4, 0.3)
                #g1_toys_pm1.SetFillStyle(3244)
                #g1_toys_pm1.SetLineColor(4)
                #g1_toys_pm1.SetLineWidth(2)
                g1_toys_pm1.SetLineWidth(0)
                g1_toys_pm1.SetOption("3")
                g1_toys_pm1.GetHistogram().SetOption("3")
                g1_toys_pm1.GetHistogram().SetOption("3")

            g1_asymptotic.SetTitle(f"Asymptotic ({limtype_label})")
            g1_asymptotic.SetFillStyle(0)
            g1_asymptotic.SetLineColor(2)
            g1_asymptotic.SetLineWidth(2)
            g1_asymptotic.SetMarkerStyle(24)
            g1_asymptotic.SetMarkerColor(2)
            g1_asymptotic.SetMarkerSize(2)
            g1_asymptotic.SetOption("LP")
            g1_asymptotic.GetHistogram().SetOption("LP")
            
            g1_toys.SetTitle(f"Toys ({limtype_label})")
            g1_toys.SetFillStyle(0)
            g1_toys.SetLineColor(4)
            g1_toys.SetLineWidth(2)
            g1_toys.SetMarkerStyle(26)
            g1_toys.SetMarkerColor(4)
            g1_toys.SetMarkerSize(2)
            g1_toys.SetOption("LP")
            g1_toys.GetHistogram().SetOption("LP")

            g1_asymptotic.GetXaxis().SetLimits(xmin, xmax)
            g1_toys.GetXaxis().SetLimits(xmin, xmax)
            
            #h1_asymptotic = cmut.root_TGraph_to_TH1(g1_asymptotic, evalBins = False, setError = False)
            #h1_toys = cmut.root_TGraph_to_TH1(g1_toys, evalBins = False, setError = False)
            
            h1_asymptotic = cmut.root_TGraph_to_TH1(g1_asymptotic, evalBins = True, setError = False)
            h1_toys = cmut.root_TGraph_to_TH1(g1_toys, evalBins = True, setError = False)
            
            h1_asymptotic.SetDirectory(0)
            h1_toys.SetDirectory(0)
            
            #h1_asymptotic.SetTitle(f"Asymptotic ({limtype})")
            h1_asymptotic.SetLineWidth(2)
            h1_asymptotic.SetLineColor(2)
            h1_asymptotic.SetMarkerSize(0)
            h1_asymptotic.SetMarkerColor(2)
            h1_asymptotic.SetFillStyle(0)
            h1_asymptotic.SetOption("hist")
            
            #h1_toys.SetTitle(f"Toys ({limtype})")
            h1_toys.SetLineWidth(2)
            h1_toys.SetLineColor(4)
            h1_toys.SetMarkerSize(0)
            h1_toys.SetMarkerColor(4)
            h1_toys.SetFillStyle(0)
            h1_toys.SetOption("hist")
            
            ymin = min(min(numpy.array(g1_asymptotic.GetY())), min(numpy.array(g1_toys.GetY())))
            ymin = 10**(round(numpy.log10(ymin))-1)
            
            ymax = max(max(numpy.array(g1_asymptotic.GetY())), max(numpy.array(g1_toys.GetY())))
            ymax = 10**(round(numpy.log10(ymax))+3)
            
            ymax_ratio = 2.0
            h1_tmp = h1_toys.Clone(f"{h1_toys.GetName()}_tmp")
            h1_tmp.SetDirectory(0)
            h1_tmp.Divide(h1_asymptotic)
            # For some crazy weird reason, the plotted histogram sometimes shows completely wrong values if this is not done.
            # However, hist.Print("all") prints the correct values
            #h1_tmp.Scale(1.0)
            #h1_tmp.Print("all")
            ymax_ratio = h1_tmp.GetBinContent(h1_tmp.GetMaximumBin())
            print(f"ymax_ratio = {ymax_ratio}")
            ymax_ratio = max(2, numpy.ceil(ymax_ratio))
            print(f"ymax_ratio = {ymax_ratio}")
            
            outfname = f"{outdir}/{args.toystamp}/compare_limits_{limtype}_{ctau}.pdf"
            
            ctau_num = ctau.replace("p", ".").split("mm")[0]
            
            h1_dummy = g1_toys.GetHistogram().Clone(f"h1_dummy_{ctau}")
            h1_dummy.SetTitle("")
            h1_dummy.Reset()
            
            utils.root_plot1D_legacy(
                #l_hist = [h1_toys, h1_asymptotic],
                l_hist = [h1_dummy],
                l_graph_overlay = [g1_toys, g1_asymptotic] +
                    [g1_toys_pm1]*int(g1_toys_pm1 is not None) +
                    [g1_asymptotic_pm1]*int(g1_asymptotic_pm1 is not None),
                ratio_num_den_pairs = [[h1_toys, h1_asymptotic]],
                #ratio_num_den_pairs = [(h1_asymptotic, h1_toys)],
                #l_ratio_hist_overlay = [h1_tmp],
                outfile = outfname,
                xrange = (xmin, xmax),
                yrange = (ymin, ymax),
                #yrange = (0, 2),
                logx = False,
                #logy = False,
                logy = True,
                ytitle = "Cross Section [fb]",
                xtitle_ratio = "m_{#tilde{#tau}} [GeV]",
                ytitle_ratio = "Toys / Asymp.",
                yrange_ratio = (0, ymax_ratio),
                ratiodrawopt = "L", # Using hist plots incorrect values for some histograms, for some strange reason
                centertitlex = True, centertitley = True,
                centerlabelx = False, centerlabely = False,
                gridx = True, gridy = True,
                ndivisionsx = None,
                #ndivisionsy_ratio = (int(2*ymax_ratio), 5, 0),
                ndivisionsy_ratio = (int(ymax_ratio), 5, 0),
                stackdrawopt = "nostack",
                legendpos = "UL",
                legendtitle = f"#splitline{{    }}{{#splitline{{{scenario}}}{{c#tau_{{0}} = {ctau_num} mm}}}}",
                legendncol = 2,
                legendtextsize = 0.05,
                legendwidthscale = 2,
                legendheightscale = 1.2,
                lumiText = "2018 (13 TeV)",
                CMSextraText = "Preliminary"
            )
    
    
    return 0

if __name__ == "__main__":
    main()