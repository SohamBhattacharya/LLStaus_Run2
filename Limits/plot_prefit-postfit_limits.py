#!/usr/bin/env -S python3 -u

import argparse
import numpy
import os
import ROOT

import utils.commonutils as cmut
import utils.utils as utils


def main () :
    
    # Argument parser
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        "--input",
        help = "Input postfit root files",
        type = str,
        nargs = "+",
        required = True,
    )
    
    parser.add_argument(
        "--eras",
        help = "Eras",
        type = str,
        required = True,
        nargs = "+",
        choices = ["2016_preVFP", "2016_postVFP", "2017", "2018", "added"],
    )
    
    parser.add_argument(
        "--channels",
        help = "Channels",
        type = str,
        required = True,
        nargs = "+",
        choices = ["BRT2",],
    )
    
    parser.add_argument(
        "--type",
        help = "background-only (b) or signal+background (s)",
        type = str,
        required = True,
        choices = ["s", "b"],
    )
    
    parser.add_argument(
        "--nosig",
        help = "Will not plot signal in the plots",
        action = "store_true",
        default = False,
    )
    
    parser.add_argument(
        "--outdir",
        help = "Output directory; if not provided, will output to the directory of the input root file",
        type = str,
        required = False,
    )
    
    parser.add_argument(
        "--title",
        help = "Title",
        type = str,
        required = False,
        default = "",
    )
    
    parser.add_argument(
        "--cmsextratext",
        help = "CMS extra text",
        type = str,
        required = False,
        default = "Private Work",
        nargs = "?",
    )
    
    
    # Parse arguments
    args = parser.parse_args()
    
    d_lumitexts = {
        "2016_preVFP": "19.5 fb^{-1} (13 TeV)",
        "2016_postVFP": "16.8 fb^{-1} (13 TeV)",
        "2017": "41.5 fb^{-1} (13 TeV)",
        "2018": "59.8 fb^{-1} (13 TeV)",
        "added": "138 fb^{-1} (13 TeV)",
    }
    
    d_channel_info = {}
    
    d_channel_info["BRT2"] = {
        "label": "BRT2",
        "nbins": 8,
        "binstart": 3,
        "eras": args.eras,
        "data": "data_obs",
        "processes_bkg": ["misid"],
        "processes_sig": ["sig"],
        "bkg_total": "TotalBkg",
        "sig_total": "TotalSig",
        "procs_total": "TotalProcs",
        "fits": ["prefit", "postfit"],
    }
    
    d_proc_info = {
        "sig": {
            "label": "",
            "color": "#bd1f01",
        },
        
        "misid": {
            "label": "Background",
            "color": "#92dadd",
        },
    }
    
    for infile_name in args.input :
        
        outdir = args.outdir if (args.outdir) else os.path.dirname(infile_name)
        os.system(f"mkdir -p {outdir}")
        
        parsed_result = cmut.parse_string_regex(
            s = infile_name,
            regexp = "SMS-TStauStau_MStau-(?P<mstau>\d+)_ctau-(?P<ctau>\d+)mm_mLSP-(?P<mlsp>\d+)",
        )
        
        mstau_str = parsed_result["mstau"]
        ctau_str = parsed_result["ctau"]
        
        #mstau = float(mstau_str)
        #ctau = float(ctau_str.replace("p", "."))
        
        sig_label = f"m_{{#tilde{{#tau}}}}={mstau_str} GeV, c#tau_{{0}}={ctau_str} mm"
        d_proc_info["sig"]["label"] = sig_label
        
        print(f"Opening file: {infile_name}")
        infile = ROOT.TFile.Open(infile_name)
        
        for ch_name in args.channels :
            
            ch_info = d_channel_info[ch_name]
            
            for fit in ch_info["fits"] :
                
                for era in ch_info["eras"] :
                    
                    d_hist = {}
                    
                    procs = ch_info["processes_bkg"] + ch_info["processes_sig"] + [
                        ch_info["data"],
                        ch_info["bkg_total"],
                        ch_info["sig_total"],
                        ch_info["procs_total"],
                    ]
                    
                    for iproc, proc in enumerate(procs) :
                        #bin9_BRT2_2016_preVFP
                        nbins = ch_info["nbins"]
                        binstart = ch_info["binstart"]
                        binend = binstart + nbins
                        hist_proc_name = f"{ch_name}_{era}_{fit}_{proc}"
                        hist_proc_name = hist_proc_name.replace("__", "_")
                        hist_proc = ROOT.TH1F(hist_proc_name, proc, nbins, 1, nbins+1)
                        hist_proc.SetDirectory(0)
                        
                        #color = utils.ColorIterator(iproc, 0)
                        color = 1
                        label = hist_proc.GetTitle()
                        
                        if (proc in d_proc_info) :
                            
                            color = ROOT.TColor.GetColor(d_proc_info[proc]["color"])
                            label = d_proc_info[proc]["label"]
                        
                        for ibin in range(nbins) :
                            binnum = ibin + binstart
                            inhist_name = f"bin{binnum}_{ch_name}_{era}_{fit}/{proc}"
                            inhist_name = inhist_name.replace("__", "_")
                            #print(inhist_name)
                            print(f"Getting histogram: {inhist_name}")
                            inhist = infile.Get(inhist_name).Clone()
                            inhist.SetDirectory(0)
                            
                            bin_val = inhist.GetBinContent(1)
                            bin_err = inhist.GetBinError(1)
                            
                            if numpy.isnan(bin_val) :
                                bin_val = 0
                                bin_err = 0
                            
                            hist_proc.SetBinContent(ibin+1, bin_val)
                            hist_proc.SetBinError(ibin+1, bin_err)
                            
                        
                        hist_proc.SetTitle(label)
                        hist_proc.SetOption("hist")
                        hist_proc.SetFillColor(color)
                        hist_proc.SetFillStyle(1001)
                        hist_proc.SetLineColor(color)
                        hist_proc.SetMarkerColor(color)
                        hist_proc.SetLineWidth(0)
                        #hist_proc.SetMarkerSize(0)
                        
                        if ("sig" in proc) :
                            
                            hist_proc.SetMarkerSize(0)
                            hist_proc.SetFillStyle(0)
                            hist_proc.SetFillColor(0)
                            hist_proc.SetLineWidth(2)
                            hist_proc.SetLineStyle(7)
                            #hist_proc.SetLineColor(color)
                            #hist_proc.SetMarkerColor(color)
                            #hist_proc.SetOption("hist")
                            
                            sig_with_bkg = True
                        
                        print(hist_proc_name)
                        hist_proc.Print()
                        
                        d_hist[proc] = hist_proc#.Clone()
                    
                    logy = True
                    ymin = 1e-2
                    hist_data = d_hist[ch_info["data"]]
                    gr_data = ROOT.TGraphAsymmErrors(hist_data.GetNbinsX())
                    
                    for bin in range(hist_data.GetNbinsX()) :
                        
                        bin_center = hist_data.GetBinCenter(bin+1)
                        bin_val = hist_data.GetBinContent(bin+1)
                        
                        if not bin_val :
                            
                            # ROOT sets the ymin to a lower value thean what is provided, when using log scale 
                            # Set the data point at the edge of ymin
                            if logy :
                                
                                hist_data.SetBinContent(bin+1, 2.51*ymin/10)
                            
                            hist_data.SetBinError(bin+1, -numpy.log((1-0.68)/2))
                        
                        unc_d, unc_u = cmut.get_garwood_errors(n = bin_val)
                        
                        gr_data.SetPoint(bin, bin_center, bin_val)
                        gr_data.SetPointError(
                            bin,
                            0.0, 0.0,
                            unc_d, unc_u
                        )
                    
                    l_hist_bkg = [d_hist[_proc] for _proc in ch_info["processes_bkg"]]
                    l_hist_sig = [d_hist[_proc] for _proc in ch_info["processes_sig"]]
                    
                    l_hist = l_hist_bkg
                    hist_bkg_total = d_hist[ch_info["bkg_total"]].Clone()
                    
                    #if (sig_with_bkg) :
                    #    
                    #    l_hist = l_hist_bkg + l_hist_sig
                    #    
                    #    hist_bkg_total = d_hist[ch_info["procs_total"]]
                    
                    hist_bkg_total.SetFillColorAlpha(1, 1.0)
                    hist_bkg_total.SetFillStyle(3354)
                    hist_bkg_total.SetMarkerSize(0)
                    hist_bkg_total.SetMarkerColor(1)
                    #hist_bkg_total.SetMarkerStyle(0)
                    hist_bkg_total.SetLineWidth(0)
                    #hist_bkg_total.SetLineColor(2)
                    hist_bkg_total.SetOption("E2")
                    
                    hist_bkg_total.SetTitle("Background uncertainty")
                    
                    gr_bkg_error = ROOT.TGraphAsymmErrors(hist_bkg_total.GetNbinsX())
                    
                    for ibin in range(hist_bkg_total.GetNbinsX()) :
                        
                        bin_val = hist_bkg_total.GetBinContent(ibin+1)
                        bin_err = hist_bkg_total.GetBinError(ibin+1)
                        bin_halfwidth = hist_bkg_total.GetBinWidth(ibin+1) / 2.0
                        
                        gr_bkg_error.SetPoint(ibin, hist_bkg_total.GetBinCenter(ibin+1), bin_val)
                        
                        gr_bkg_error.SetPointError(
                            ibin,
                            bin_halfwidth, bin_halfwidth,
                            bin_err, bin_err
                        )
                    
                    gr_bkg_error.SetName("gr_bkg_error")
                    gr_bkg_error.SetTitle("Background uncertainty")
                    gr_bkg_error.SetFillStyle(3254)
                    gr_bkg_error.SetFillColorAlpha(1, 0.5)
                    gr_bkg_error.SetMarkerColor(1)
                    gr_bkg_error.SetLineWidth(0)
                    gr_bkg_error.GetHistogram().SetOption("E2")
                    
                    bin_max = max(hist_data.GetBinContent(hist_data.GetMaximumBin()), hist_bkg_total.GetBinContent(hist_bkg_total.GetMaximumBin()))
                    
                    yrange_power = int(numpy.log10(bin_max))+1
                    
                    title = f"{fit.title()}"
                    out_suffix = ""
                    
                    if fit == "postfit" :
                        
                        out_suffix = f"_{args.type}"
                        
                        if args.type == "s" :
                            
                            title = f"{title} (sig.+bkg. hypothesis)"
                        
                        elif args.type == "b" :
                            
                            title = f"{title} (background-only hypothesis)"
                    
                    plotfile = f"{outdir}/SMS-TStauStau_MStau-{mstau_str}_ctau-{ctau_str}mm_mLSP-1_{ch_name}_{era}_{fit}{out_suffix}.pdf".replace("__", "_")
                    
                    xmin = hist_data.GetXaxis().GetXmin()
                    xmax = hist_data.GetXaxis().GetXmax()
                    
                    gr_unity = ROOT.TGraph()
                    gr_unity.AddPoint(xmin, 1)
                    gr_unity.AddPoint(xmax, 1)
                    gr_unity.SetLineColor(1)
                    gr_unity.SetLineWidth(1)
                    gr_unity.GetHistogram().SetOption("L")
                    
                    hist_data_ratio_bkg = hist_data.Clone()
                    
                    hist_data_ratio_tot = hist_data.Clone()
                    hist_data_ratio_tot.SetLineColor(2)
                    hist_data_ratio_tot.SetLineStyle(7)
                    hist_data_ratio_tot.SetMarkerColor(2)
                    hist_data_ratio_tot.SetOption("NOERR")
                    
                    ratio_num_den_pairs = [
                        (hist_data_ratio_bkg, d_hist[ch_info["bkg_total"]]),
                        (hist_data_ratio_tot, d_hist[ch_info["procs_total"]])
                    ]
                    
                    gr_ratio_num_unc = ROOT.TGraphAsymmErrors(hist_data.GetNbinsX())
                    
                    for ibin in range(hist_data.GetNbinsX()) :
                        
                        bin_center = gr_data.GetPointX(ibin)
                        num = gr_data.GetPointY(ibin)
                        num_unc_u = gr_data.GetErrorYhigh(ibin)
                        num_unc_d = gr_data.GetErrorYlow(ibin)
                        
                        den = d_hist[ch_info["bkg_total"]].GetBinContent(ibin+1)
                        den_unc_u = d_hist[ch_info["bkg_total"]].GetBinError(ibin+1)
                        den_unc_d = d_hist[ch_info["bkg_total"]].GetBinError(ibin+1)
                        
                        ratio = num/den if den else 0
                        gr_ratio_unc_u = 0
                        gr_ratio_unc_d = 0
                        
                        if num :
                            gr_ratio_unc_u = ratio*num_unc_u/num
                            gr_ratio_unc_d = ratio*num_unc_d/num
                        
                        elif not num and den :
                            gr_ratio_unc_u =  num_unc_u / den
                            gr_ratio_unc_d =  0
                        
                        gr_ratio_num_unc.SetPointX(ibin, bin_center)
                        gr_ratio_num_unc.SetPointY(ibin, ratio)
                        gr_ratio_num_unc.SetPointEYhigh(ibin, gr_ratio_unc_u)
                        gr_ratio_num_unc.SetPointEYlow(ibin, gr_ratio_unc_d)
                    
                    gr_ratio_num_unc.SetLineColor(1)
                    gr_ratio_num_unc.SetLineWidth(2)
                    gr_ratio_num_unc.SetMarkerColor(1)
                    gr_ratio_num_unc.SetMarkerSize(2)
                    gr_ratio_num_unc.SetMarkerStyle(20)
                    gr_ratio_num_unc.SetFillStyle(0)
                    gr_ratio_num_unc.GetHistogram().SetOption("PE1")
                    
                    gr_ratio_num_unc.Print("all")
                    
                    ytitle_ratio = "Data / Pred."
                    legendncol = 2
                    legendwidthscale = 1.9
                    
                    if args.nosig :
                        
                        l_hist_sig = []
                        ratio_num_den_pairs = ratio_num_den_pairs[0: 1]
                        ytitle_ratio = "Data / Bkg."
                        legendncol = 1
                        legendwidthscale = 0.9
                    
                    outfilename_root = plotfile.replace(".pdf", ".root")
                    outfile_root = ROOT.TFile.Open(outfilename_root, "RECREATE")
                    outfile_root.cd()
                    
                    for hist in [*l_hist, *l_hist_sig, hist_data] :
                        
                        hist_tmp = hist.Clone()
                        hist_tmp.SetLineWidth(2)
                        hist_tmp.SetLineColor(1)
                        hist_tmp.SetMarkerSize(2)
                        hist_tmp.SetMarkerColor(1)
                        hist_tmp.SetFillStyle(0)
                        hist_tmp.SetOption("hist E1")
                        hist_tmp.Write()
                    
                    hist_data.Reset()
                    hist_data.SetTitle("Data")
                    hist_data.SetLineColor(1)
                    hist_data.SetMarkerColor(1)
                    hist_data.SetMarkerSize(2)
                    hist_data.SetMarkerStyle(20)
                    hist_data.SetLineWidth(2)
                    hist_data.SetOption("PE1")
                    hist_data.SetFillStyle(0)
                    hist_data.SetFillColor(0)
                    
                    gr_data.SetLineColor(1)
                    gr_data.SetLineWidth(2)
                    gr_data.SetMarkerColor(1)
                    gr_data.SetMarkerSize(2)
                    gr_data.SetMarkerStyle(20)
                    gr_data.SetFillStyle(0)
                    gr_data.GetHistogram().SetOption("PE1")
                    
                    utils.root_plot1D_legacy(
                        l_hist = l_hist,
                        l_hist_overlay = [*l_hist_sig, hist_data],
                        l_graph_overlay = [gr_bkg_error, gr_data],
                        l_ratio_graph_overlay = [gr_unity, gr_ratio_num_unc],
                        l_legend_order = [hist_data, *l_hist, gr_bkg_error, *l_hist_sig],
                        outfile = plotfile,
                        xrange = (xmin, xmax),
                        yrange = (ymin, 10**(yrange_power+3)),
                        ratio_num_den_pairs = ratio_num_den_pairs,
                        ratio_mode = "data",
                        no_xerror = True,
                        logx = False,
                        logy = logy,
                        xtitle = "",
                        ytitle = "Events",
                        xtitle_ratio = "SR bins",
                        ytitle_ratio = ytitle_ratio,
                        yrange_ratio = (0, 5),
                        centertitlex = True,
                        centertitley = True,
                        centerlabelx = True,
                        centerlabely = False,
                        gridx = False,
                        gridy = False,
                        ndivisionsx = (ch_info["nbins"], 1, 0),
                        ndivisionsy = None,
                        ndivisionsy_ratio = (5, 5, 0),
                        stackdrawopt = "",
                        ratiodrawopt = "E1",
                        legendpos = "UL",
                        legendncol = legendncol,
                        legendtextsize = 0.05,
                        legendtitle = title,
                        legendheightscale = 1.2,
                        legendwidthscale = legendwidthscale,
                        legendpadleft_extra = -0.03,
                        CMSextraText = args.cmsextratext,
                        lumiText = f"{d_lumitexts[era]}"
                    )
                    
                    outfile_root.Close()
                    
                    #os.system(f"pdftoppm -png -r 600 -cropbox -singlefile {plotfile} {plotfile}")
        
        infile.Close()
    
    return 0


if __name__ == "__main__" :
    
    main()