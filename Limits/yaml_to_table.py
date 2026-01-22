#!/usr/bin/env python3

import argparse
from decimal import FloatOperation
import itertools
import functools
from matplotlib.font_manager import FontProperties
import numpy
import operator
import os
import pandas

from prettytable import PrettyTable

import ROOT
from scipy.fftpack import cc_diff
ROOT.gROOT.SetBatch(True)

import utils.utils as utils
import utils.commonutils as cmut
import utils.pdgRounding as pdgrnd


DELTA_POISSON_0 = -numpy.log((1-0.68)/2)

def main() :
    
    # Argument parser
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        "--inputs",
        help = (
            "Input yaml filenames. "
            "Optionally a key to be read from the file can be passed as input_filename:key. "
            "An output will be created for each input."
        ),
        type = str,
        nargs = "*",
        required = True,
    )
    
    parser.add_argument(
        "--outdir",
        help = "Output directory",
        type = str,
        required = True,
    )
    
    parser.add_argument(
        "--type",
        help = "Type of information in the input file",
        type = str,
        required = True,
        choices = ["yields", "systematics"],
    )
    
    parser.add_argument(
        "--config",
        help = "Configuration yaml for parsing and formatting the input",
        type = str,
        required = True,
    )
    
    parser.add_argument(
        "--tex",
        help = "Latex template file. Will use this to create a pdf and of the table. Pass empty string \"\" to disable latex compilation",
        type = str,
        required = False,
        default = "utils/tex/template_standalone_png.tex",
    )
    
    parser.add_argument(
        "--plot",
        help = "Create a plot of the yields",
        action = "store_true",
        default = False,
    )
    
    parser.add_argument(
        "--cmsextratext",
        help = "CMS extra text",
        type = str,
        required = False,
        default = "Preliminary",
        nargs = "?",
    )
    
    parser.add_argument(
        "--lumitext",
        help = "CMS lumitext",
        type = str,
        required = False,
        default = "(13 TeV)",
        nargs = "?",
    )
    
    parser.add_argument(
        "--blind",
        help = "Will not include observation; will append \"_blinded\" to the output directory",
        action = "store_true",
        default = False,
    )
    
    
    # Parse arguments
    args = parser.parse_args()
    
    if not args.cmsextratext :
        
        args.cmsextratext = ""
    
    if not args.lumitext :

        args.lumitext = ""

    d_config = cmut.load_config(args.config)
    
    for infilename in args.inputs :
        
        input_key = None
        
        if (":" in infilename) :
            
            infilename, input_key = infilename.split(":")
        
        d_input = cmut.load_config(infilename)
        d_input = d_input.get(input_key, d_input)
        print(d_input)
        
        if (args.type == "yields") :
            
            bin_header = d_config["bin_header"]
            d_table = {}
            d_table[bin_header] = []
            
            l_bins_edges = numpy.array(d_input["bin_edges"], dtype = float)
            binstart = d_input["binstart"]
            nbins = len(l_bins_edges)-1
            l_plot_bins = list(range(1, nbins+1))
            bin_args = [nbins, l_bins_edges]
            
            if ("plot_bins" in d_config) :
                
                l_plot_bins = d_config["plot_bins"]
                nbins = len(l_plot_bins)
                bin_args = [nbins, 1, nbins+1]
            
            
            l_proc_types = [] if args.blind else ["obs"]
            l_proc_types.extend(["bkg", "sig"])
            
            d_hist_yields = {}
            d_gr_uncs = {}
            
            for proc_type in l_proc_types :
                
                if (proc_type not in d_config) :
                    
                    cmut.logger.warning(f"[process type {proc_type}] not found in config. Skipping.")
                    continue
                
                if (proc_type not in d_hist_yields) :
                    
                    d_hist_yields[proc_type] = {}
                    d_gr_uncs[proc_type] = {}
                
                for proc_name, d_proc_cfg in d_config[proc_type].items() :
                    
                    if proc_name not in d_input :
                        
                        cmut.logger.warning(f"[process {proc_name}] not found in input. Skipping.")
                        continue
                    
                    proc_label = d_proc_cfg["label"]
                    d_table[proc_label] = []
                    
                    #to_plot = (proc_type in ["obs", "sig"]) or (proc_name in d_config["plot_bkgs"])
                    #to_plot = True
                    #
                    #if to_plot :
                    
                    #d_hist_yields[proc_type][proc_name] = ROOT.TH1F(f"h1_{proc_name}", d_proc_cfg.get("rootlabel", proc_label), nbins, 1, 1+nbins)
                    d_hist_yields[proc_type][proc_name] = ROOT.TH1F(f"h1_{proc_name}", d_proc_cfg.get("rootlabel", proc_label), *bin_args)
                    color = d_proc_cfg.get("color", "#00000")
                    color = ROOT.TColor.GetColor(color) if isinstance(color, str) else color
                    
                    d_hist_yields[proc_type][proc_name].SetLineColor(color)
                    d_hist_yields[proc_type][proc_name].SetLineStyle(d_proc_cfg.get("linestyle", 1))
                    d_hist_yields[proc_type][proc_name].SetMarkerColor(color)
                    
                    d_gr_uncs[proc_type][proc_name] = ROOT.TGraphAsymmErrors(nbins)
                    d_gr_uncs[proc_type][proc_name].SetName(f"gr_{proc_name}_unc")
                    
                    if proc_type == "obs" :
                        d_gr_uncs[proc_type][proc_name].SetTitle(proc_label)
                    elif proc_type == "bkg":
                        d_gr_uncs[proc_type][proc_name].SetTitle(d_config.get("unc_label_bkg", "Background uncertainty"))
                    
                    d_proc_info = d_input[proc_name]
                    
                    if "bins" not in d_config :
                        
                        d_config["bins"] = {}
                        
                        for ibin, bin in enumerate(l_plot_bins) :
                            
                            d_config["bins"][bin] = {
                                "label": f"{l_bins_edges[ibin]}--{l_bins_edges[ibin+1]}",
                            }
                        #
                        #d_config["bins"][l_plot_bins[-1]] = {
                        #    "label": f"${{>}}{l_bins_edges[-1]}$",
                        #}
                    
                    for bin, d_bin_cfg in d_config["bins"].items() :
                        
                        bin_label = d_bin_cfg["label"]
                        
                        if bin_label not in d_table[bin_header] :
                            
                            d_table[bin_header].append(bin_label)
                        
                        d_yield_info = d_proc_info[bin]
                        
                        val = d_yield_info["yield"]
                        
                        if (proc_type == "obs") :
                            
                            d_table[proc_label].append(f"{int(val)}")
                            unc = val**0.5
                            unc_d, unc_u = cmut.get_garwood_errors(n = val)
                            
                            # Prediction can be a treated as fake observation 
                            if "unc_stat" in d_yield_info and "unc_syst" in d_yield_info :
                                
                                unc_stat = d_yield_info["unc_stat"]
                                unc_syst = d_yield_info["unc_syst"]
                                unc = (unc_stat**2 + unc_syst**2)**0.5
                                unc_d = unc
                                unc_u = unc
                        
                        else :
                            
                            unc_stat = d_yield_info["unc_stat"]
                            unc_syst = d_yield_info["unc_syst"]
                            unc = (unc_stat**2 + unc_syst**2)**0.5
                            
                            if val :
                                
                                #val_str, _ = pdgrnd.pdgRound(val, 0)
                                #_, unc_stat_str = pdgrnd.pdgRound(val, unc_stat)# if (unc_stat > 0.05) else (None, f"{unc_stat: 0.1f}")
                                #_, unc_syst_str = pdgrnd.pdgRound(val, unc_syst)# if (unc_syst > 0.05) else (None, f"{unc_syst: 0.1f}")
                                
                                #val_str = f"{val:0.2g}"# if val < 99.5 else f"{val:0.0f}"
                                #unc_stat_str = f"{unc_stat:0.2g}"# if unc_stat < 99.5 else f"{unc_stat:0.0f}"
                                #unc_syst_str = f"{unc_syst:0.2g}"# if unc_syst < 99.5 else f"{unc_syst:0.0f}"
                                
                                val_str = f"{val:0.4}"# if val < 99.5 else f"{val:0.0f}"
                                unc_stat_str = f"{unc_stat:0.4g}"# if unc_stat < 99.5 else f"{unc_stat:0.0f}"
                                unc_syst_str = f"{unc_syst:0.4g}"# if unc_syst < 99.5 else f"{unc_syst:0.0f}"
                                
                                #val_str, unc_stat_str, unc_syst_str = cmut.cms_rounding(val, unc_stat, unc_syst)
                                
                                d_table[proc_label].append(f"${val_str}{{\\pm}}{unc_stat_str}{{\\pm}}{unc_syst_str}$")
                            
                            else :
                                
                                _, unc_str = pdgrnd.pdgRound(val, unc)# if (unc > 0.05) else (None, f"{unc: 0.1f}")
                                
                                d_table[proc_label].append(f"${{<}}{unc_str}$")
                        
                        # Histograms and graphs
                        #if to_plot and bin in l_plot_bins :
                        if bin in l_plot_bins :
                            
                            bin_num = l_plot_bins.index(bin)
                            d_hist_yields[proc_type][proc_name].SetBinContent(bin_num+1, val)
                            d_hist_yields[proc_type][proc_name].SetBinError(bin_num+1, unc)
                            
                            bin_center = d_hist_yields[proc_type][proc_name].GetXaxis().GetBinCenter(bin_num+1)
                            bin_halfwidth = 0.5*d_hist_yields[proc_type][proc_name].GetXaxis().GetBinWidth(bin_num+1)
                            
                            if (proc_type == "obs" and not d_proc_cfg.get("fakeobs", False)) :
                                
                                d_gr_uncs[proc_type][proc_name].SetPoint(bin_num, bin_center, val)
                                
                                d_gr_uncs[proc_type][proc_name].SetPointError(
                                    bin_num,
                                    0.0, 0.0,
                                    unc_d, unc_u
                                )
                            
                            else :
                                
                                d_gr_uncs[proc_type][proc_name].SetPoint(bin_num, bin_center, val)
                                
                                d_gr_uncs[proc_type][proc_name].SetPointError(
                                    bin_num,
                                    bin_halfwidth, bin_halfwidth,
                                    (unc if val else 0), unc
                                )
            
            if len(d_config["bkg"]) <= 2 and "total_bkg" in d_config["bkg"] :
                
                total_bkg_label = d_config["bkg"]["total_bkg"]["label"]
                del d_table[total_bkg_label]
        
        elif (args.type == "systematics") :
            
            syst_header = d_config["syst_header"]
            d_table = {}
            d_table[syst_header] = []
            
            for proc_name, d_proc_cfg in d_config["procs"].items() :
                
                proc_label = d_proc_cfg["label"]
                d_table[proc_label] = []
                
                d_proc_info = d_input[proc_name]
                
                for syst_name, d_syst_cfg in d_config["systematics"].items() :
                    
                    syst_label = d_syst_cfg["label"]
                    
                    if syst_label not in d_table[syst_header] :
                        
                        d_table[syst_header].append(syst_label)
                    
                    # If the process does not have this uncertainty, add a dummy entry
                    if syst_name not in d_proc_info :
                        
                        d_table[proc_label].append("--")
                        continue
                    
                    d_syst_info = d_proc_info[syst_name]
                    l_syst_vals = []
                    
                    for bin, syst_info in d_syst_info.items() :
                        
                        if (syst_info["unc_rel"] > 0) :
                            
                            l_syst_vals.append(syst_info["unc_rel"])
                    
                    if not l_syst_vals :
                        
                        cmut.logger.warning(f"[process {proc_name}] [systematic {syst_name}] has no valid entries.")
                        
                        d_table[proc_label].append("--")
                        continue
                    
                    syst_min = 100*min(l_syst_vals)
                    syst_max = 100*max(l_syst_vals)
                    syst_med = 100*numpy.median(l_syst_vals)
                    
                    #syst_min_str, _ = pdgrnd.pdgRound(syst_min, 0)
                    #syst_max_str, _ = pdgrnd.pdgRound(syst_max, 0)
                    #syst_med_str, _ = pdgrnd.pdgRound(syst_med, 0)
                    
                    syst_min_str = f"{syst_min:0.2g}" if syst_min < 99.5 else f"{syst_min:0.0f}"
                    syst_max_str = f"{syst_max:0.2g}" if syst_max < 99.5 else f"{syst_max:0.0f}"
                    syst_med_str = f"{syst_med:0.2g}" if syst_med < 99.5 else f"{syst_med:0.0f}"
                    
                    #syst_med_str, syst_min_str, syst_max_str = cmut.cms_rounding(syst_med, syst_min, syst_max)
                    
                    syst_min_str = str(syst_min_str) if (float(syst_min_str) >= 0.1) else "{<}0.1"
                    syst_max_str = str(syst_max_str) if (float(syst_max_str) >= 0.1) else "{<}0.1"
                    syst_med_str = str(syst_med_str) if (float(syst_med_str) >= 0.1) else "{<}0.1"
                    
                    # Set to the median if they differ by less than 1
                    #if abs(syst_min-syst_med) < 0.5 and abs(syst_max-syst_med) < 0.5 :
                    if abs(syst_max-syst_min) <= 1 :
                        
                        d_table[proc_label].append(f"${syst_med_str}\\%$")
                    
                    else :
                        
                        d_table[proc_label].append(f"${syst_min_str},\\ {syst_med_str},\\ {syst_max_str}\\%$")
            
            #print(d_table)
        
        
        # Print table to console
        ptable = PrettyTable()
        
        for header, column in d_table.items() :
            
            l_replace_pairs = [
                ("{\\pm}", " ± "),
                ("$", ""),
                ("{<}", "<"),
                ("\\%", "%"),
                ("\\ ", " "),
                ("{", ""),
                ("}", ""),
                ("[", ""),
                ("]", ""),
                ("\\tau", "τ"),
                ("\\PSGt", "stau"),
            ]
            
            for replace_pair in l_replace_pairs :
                
                header = header.replace(*replace_pair)
                column = [str(_element).replace(*replace_pair) for _element in column]
                
            ptable.add_column(header, column)
        
        print("\n")
        print(ptable)
        
        
        # Print and save latex table
        latex_table_str = pandas.DataFrame.from_dict(d_table, orient = "columns", dtype = str).to_latex(
            index = False,
            column_format = "l"+"c"*(len(d_table)-1)
        )
        
        l_replace_pairs = [
            #("toprule", "hline"),
            ("\\toprule", ""),
            ("midrule", "hline"),
            #("\\midrule", ""),
            ("midrule", "hline"),
            ("bottomrule", "hline"),
            ("\\begin{tabular}", "\\fittable{\\begin{tabular}"),
            ("\\end{tabular}", "\\end{tabular}}"),
        ]
        
        for replace_pair in l_replace_pairs :
            
            latex_table_str = latex_table_str.replace(*replace_pair)
        
        print("\n")
        print(latex_table_str)
        
        #outdir = args.outdir if args.outdir else os.path.dirname(infilename)
        outdir = f"{args.outdir}_blinded" if args.blind else args.outdir
        os.system(f"mkdir -p {outdir}")
        
        infilename_noext = os.path.splitext(os.path.basename(infilename))[0]
        
        outfilename = f"{outdir}/{infilename_noext}.tex"
        cmut.logger.info(f"Writing table to: {outfilename}")
        
        with open(outfilename, "w") as fopen :
            
            fopen.write(latex_table_str)
        
        
        if (args.tex) :
            
            # Standalone png and pdf
            tex_template = ""
            latex_table_full_str = f"\\begin{{table}}\n\n{latex_table_str}\n\n\\end{{table}}"
            
            with open(args.tex, "r") as fopen :
                
                tex_template = fopen.read()
            
            #print(tex_template)
            tex_template = tex_template.replace("%CONTENT%", latex_table_full_str)
            
            outfilename = f"{outdir}/{infilename_noext}_standalone.tex"
            outfilename_noext, _ = os.path.splitext(outfilename)
            cmut.logger.info(f"Writing table to: {outfilename}")
            
            with open(outfilename, "w") as fopen :
                
                fopen.write(tex_template)
            
            retval = os.system(f"latexmk -shell-escape -silent -jobname={outfilename_noext} {outfilename}; cd {outdir}; latexmk -c")
            
            if retval :
                
                cmut.logger.error("Error compiling latex.")
                return retval
        
        
        if (args.plot) :
            
            has_obs = "obs" in d_hist_yields
            obs_key = d_config.get("obs_key", "obs")
            hist_obs = d_hist_yields["obs"][obs_key] if has_obs else None
            gr_obs_unc = d_gr_uncs["obs"][obs_key] if has_obs else None
            total_bkg_key = d_config.get("total_bkg_key", "total_bkg")
            
            if has_obs :
                
                hist_ratio = d_hist_yields["obs"][obs_key].Clone("ratio")
                gr_ratio_num_unc = d_gr_uncs["obs"][obs_key].Clone("ratio_num_unc")
            
            else :
                
                hist_ratio = d_hist_yields["bkg"][total_bkg_key].Clone("ratio")
                gr_ratio_num_unc = d_gr_uncs["bkg"][total_bkg_key].Clone("ratio_num_unc")
            
            gr_bkg_unc = d_gr_uncs["bkg"][total_bkg_key]
            gr_ratio_unc = d_gr_uncs["bkg"][total_bkg_key].Clone("ratio_unc")
            
            for bin in range(hist_ratio.GetNbinsX()) :
                
                num = hist_ratio.GetBinContent(bin+1)
                #num_unc = hist_ratio.GetBinError(bin+1)
                num_unc_u = gr_ratio_num_unc.GetErrorYhigh(bin)
                num_unc_d = gr_ratio_num_unc.GetErrorYlow(bin)
                
                den = d_hist_yields["bkg"][total_bkg_key].GetBinContent(bin+1)
                den_unc_u = d_gr_uncs["bkg"][total_bkg_key].GetErrorYhigh(bin)
                den_unc_d = d_gr_uncs["bkg"][total_bkg_key].GetErrorYlow(bin)
                
                ratio = num/den if den else 0
                gr_ratio_unc_u = 0
                gr_ratio_unc_d = 0
                ratio_unc = 0
                
                if num :
                    
                    #ratio_unc = ratio*num_unc/num
                    gr_ratio_unc_u = ratio*num_unc_u/num
                    gr_ratio_unc_d = ratio*num_unc_d/num
                
                elif not num and den :
                    
                    # One sided Poisson error
                    #ratio_unc = -numpy.log((1-0.68)/2) / den
                    gr_ratio_unc_u =  num_unc_u / den
                    gr_ratio_unc_d =  0
                
                hist_ratio.SetBinContent(bin+1, ratio)
                #hist_ratio.SetBinError(
                #    bin+1,
                #    #ratio*num_unc/num if num else 0
                #    ratio_unc
                #)
                
                gr_ratio_num_unc.SetPointY(bin, ratio)
                gr_ratio_num_unc.SetPointEYhigh(bin, gr_ratio_unc_u)
                gr_ratio_num_unc.SetPointEYlow(bin, gr_ratio_unc_d)
                
                gr_ratio_unc.SetPointY(bin, 1)
                gr_ratio_unc.SetPointEYhigh(bin, den_unc_u/den if den else 0)
                gr_ratio_unc.SetPointEYlow(bin, den_unc_d/den if den else 0)
            
            #hist_ratio.Fit(
            #    "pol1",
            #    option = "SEM",
            #    goption = "L",
            #)
            
            for hist in [hist_obs, hist_ratio] :
                
                if hist is None :
                    continue
                
                # Only use this to draw the legend marker
                # Otherwise, the error bar style of the TGraphAsymmErrors do not match up between the plot and the legend
                hist.Reset()
                
                hist.SetLineColor(1)
                hist.SetLineWidth(2)
                #hist.SetLineWidth(0)
                hist.SetMarkerColor(1)
                hist.SetMarkerSize(2)
                hist.SetMarkerStyle(20)
                hist.SetFillStyle(0)
                hist.SetOption("PE1")
                #hist.SetOption("P")
            
            for gr in [gr_obs_unc, gr_ratio_num_unc] :
                
                if gr is None :
                    continue
                
                gr.SetLineColor(1)
                gr.SetLineWidth(2)
                gr.SetMarkerColor(1)
                gr.SetMarkerSize(2)
                gr.SetMarkerStyle(20)
                gr.SetFillStyle(0)
                gr.GetHistogram().SetOption("PE1")
                #gr.SetOption("PE1")
            
            for _, hist in d_hist_yields["bkg"].items() :
                
                hist.SetLineWidth(0)
                hist.SetMarkerSize(0)
                hist.SetFillColor(hist.GetLineColor())
                hist.SetOption("hist")
            
            for _, hist in d_hist_yields["sig"].items() :
                
                hist.SetLineWidth(3)
                hist.SetMarkerSize(0)
                hist.SetOption("hist")
            
            for gr in [gr_bkg_unc, gr_ratio_unc] :
                
                gr.SetFillStyle(3254)
                #gr.SetFillStyle(3002)
                gr.SetFillColorAlpha(1, 0.5)
                gr.SetMarkerColor(1)
                gr.SetMarkerSize(0)
                gr.SetLineWidth(0)
                gr.GetHistogram().SetOption("E2")
            
            xmin = d_hist_yields["bkg"][total_bkg_key].GetXaxis().GetXmin()
            xmax = d_hist_yields["bkg"][total_bkg_key].GetXaxis().GetXmax()
            
            ymin = 1e-2
            ymax = 10**(round(numpy.log10(d_hist_yields["bkg"][total_bkg_key].GetMaximum()))+3)
            
            logy = d_config.get("logy", False)
            
            for bin in range(hist_obs.GetNbinsX()) :
                
                bin_val = hist_obs.GetBinContent(bin+1)
                
                if not bin_val :
                    
                    # ROOT sets the ymin to a lower value thean what is provided, when using log scale 
                    # Set the data point at the edge of ymin
                    if logy :
                        
                        hist_obs.SetBinContent(bin+1, 2.51*ymin/10)
                    
                    hist_obs.SetBinError(bin+1, -numpy.log((1-0.68)/2))
            
            gr_unity = ROOT.TGraph()
            gr_unity.AddPoint(xmin, 1)
            gr_unity.AddPoint(xmax, 1)
            gr_unity.SetLineColor(1)
            gr_unity.SetLineWidth(1)
            gr_unity.GetHistogram().SetOption("L")
            
            outfilename = f"{outdir}/{infilename_noext}_plot.pdf"
            
            #hist_dummy =  ROOT.TH1F(f"h1_dummy", "h1_dummy", nbins, 1, 1+nbins)
            hist_dummy =  ROOT.TH1F(f"h1_dummy", "h1_dummy", *bin_args)
            
            l_hist_obs = [hist_obs] if has_obs else []
            #l_hist_obs = []
            l_hist_bkg = [d_hist_yields["bkg"][_proc] for _proc in d_config["plot_bkgs"]]
            l_hist_sig = list(d_hist_yields["sig"].values())
            
            legend_order_bkgs = d_config.get("legend_order_bkgs", d_config["plot_bkgs"])
            legend_order_bkgs = [d_hist_yields["bkg"][_proc] for _proc in legend_order_bkgs]
            
            l_obs_bkg_legend = (l_hist_obs + legend_order_bkgs + [gr_bkg_unc]).copy()
            #l_obs_bkg_legend = ([gr_obs_unc] + legend_order_bkgs + [gr_bkg_unc]).copy()
            l_hist_sig_legend = l_hist_sig.copy()
            
            if len(l_obs_bkg_legend) and len(l_hist_sig_legend) :
                
                nextra = len(l_obs_bkg_legend) - len(l_hist_sig_legend)
                
                if nextra < 0 :
                    
                    l_obs_bkg_legend.extend([0]*abs(nextra))
                
                elif nextra > 0 :
                    
                    l_hist_sig_legend.extend([0]*nextra)
            
            #print(l_hist_sig)
            #print(l_hist_bkg)
            #print(l_hist_obs)
            
            #hist_ratio.Print("all")
            
            yrange = d_config.get("yrange", (None, None))
            yrange = (
                yrange[0] if yrange[0] is not None else ymin,
                yrange[1] if yrange[1] is not None else ymax,
            )
            
            utils.root_plot1D_legacy(
                l_hist = l_hist_bkg,
                #l_hist_overlay = l_hist_sig + l_hist_obs,
                l_hist_overlay = l_hist_sig,
                ratio_num_den_pairs = [],
                l_graph_overlay = [gr_bkg_unc] + [gr_obs_unc],
                l_ratio_hist_overlay = [hist_ratio] if has_obs else [hist_dummy],
                l_ratio_graph_overlay = [gr_unity, gr_ratio_unc, gr_ratio_num_unc],
                #l_legend_order = l_hist_obs + [0]*int(has_obs) + l_hist_bkg + [gr_bkg_unc] + l_hist_sig,
                l_legend_order = list(itertools.chain(*zip(l_obs_bkg_legend, l_hist_sig_legend))) if d_config.get("group_sigs", False) else l_obs_bkg_legend + l_hist_sig_legend,
                outfile = outfilename,
                #xrange = (xmin, xmax),
                xrange = d_config.get("xrange", (xmin, xmax)),
                yrange = yrange,
                logx = False,
                logy = logy,
                #xtitle = "Bin number",
                ytitle = d_config.get("ytitle", "Events"),
                #xtitle_ratio = "Bin number",
                xtitle_ratio = d_config.get("bin_header_root", d_config["bin_header"].strip("$")),
                ytitle_ratio = d_config.get("ytitle_ratio", "Data / Bkg."),
                yrange_ratio = d_config.get("yrange_ratio", (0, 3)),
                centertitlex = True, centertitley = True,
                centerlabelx = d_config.get("centerlabelx", False),
                centerlabely = False,
                gridx = False, gridy = False,
                ndivisionsx = d_config.get("ndivisionsx", (nbins, 1, 0)),
                stackdrawopt = "",
                #ratiodrawopt = "E1",
                ndivisionsy_ratio = d_config.get("ndivisionsy_ratio", (3, 5, 0)),
                legendpos = "UL",
                legendtitle = d_config.get("legendtitle", ""),
                legendncol = 2,
                legendtextsize = 0.05,
                legendwidthscale = 1.9,
                legendheightscale = d_config.get("legendheightscale", 0.9),
                legendpadleft_extra = -0.03,
                no_xerror = True,
                CMSextraText = args.cmsextratext,
                lumiText = args.lumitext
            )
        
        
    return 0


if (__name__ == "__main__") :
    
    main()