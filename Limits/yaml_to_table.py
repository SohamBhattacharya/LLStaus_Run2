#!/usr/bin/env python3

import argparse
from decimal import FloatOperation
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
    )
    
    parser.add_argument(
        "--lumitext",
        help = "CMS lumitext",
        type = str,
        required = False,
        default = "(13 TeV)",
    )
    
    parser.add_argument(
        "--blind",
        help = "Will not include observation; will append \"_blinded\" to the output directory",
        action = "store_true",
        default = False,
    )
    
    
    # Parse arguments
    args = parser.parse_args()
    
    d_config = cmut.load_config(args.config)
    
    for infilename in args.inputs :
        
        input_key = None
        
        if (":" in infilename) :
            
            infilename, input_key = infilename.split(":")
        
        d_input = cmut.load_config(infilename)
        d_input = d_input.get(input_key, d_input)
        #print(d_input)
        
        if (args.type == "yields") :
            
            bin_header = d_config["bin_header"]
            d_table = {}
            d_table[bin_header] = []
            
            d_hist_yields = {}
            d_gr_uncs = {}
            
            l_plot_bins = d_config["plot_bins"]
            nbins = len(l_plot_bins)
            
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
                    
                    d_hist_yields[proc_type][proc_name] = ROOT.TH1F(f"h1_{proc_name}", d_proc_cfg.get("rootlabel", proc_label), nbins, 1, 1+nbins)
                    color = ROOT.TColor.GetColor(d_proc_cfg.get("color", "#00000"))
                    
                    d_hist_yields[proc_type][proc_name].SetLineColor(color)
                    d_hist_yields[proc_type][proc_name].SetLineStyle(d_proc_cfg.get("linestyle", 1))
                    d_hist_yields[proc_type][proc_name].SetMarkerColor(color)
                    
                    d_gr_uncs[proc_type][proc_name] = ROOT.TGraphAsymmErrors(nbins)
                    d_gr_uncs[proc_type][proc_name].SetName(f"gr_{proc_name}_unc")
                    d_gr_uncs[proc_type][proc_name].SetTitle("Bkg. uncertainty")
                    
                    d_proc_info = d_input[proc_name]
                    
                    for bin, d_bin_cfg in d_config["bins"].items() :
                        
                        bin_label = d_bin_cfg["label"]
                        
                        if bin_label not in d_table[bin_header] :
                            
                            d_table[bin_header].append(bin_label)
                        
                        d_yield_info = d_proc_info[bin]
                        
                        val = d_yield_info["yield"]
                        
                        if (proc_type == "obs") :
                            
                            d_table[proc_label].append(f"{int(val)}")
                            unc = val**0.5
                        
                        else :
                            
                            unc_stat = d_yield_info["unc_stat"]
                            unc_syst = d_yield_info["unc_syst"]
                            unc = (unc_stat**2 + unc_syst**2)**0.5
                            
                            #if val >= 0.05 :
                            if val :
                                
                                val_str, _ = pdgrnd.pdgRound(val, 0)
                                _, unc_stat_str = pdgrnd.pdgRound(val, unc_stat)# if (unc_stat > 0.05) else (None, f"{unc_stat: 0.1f}")
                                _, unc_syst_str = pdgrnd.pdgRound(val, unc_syst)# if (unc_syst > 0.05) else (None, f"{unc_syst: 0.1f}")
                                
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
                            
                            if (proc_type != "obs") :
                                
                                bin_center = d_hist_yields[proc_type][proc_name].GetXaxis().GetBinCenter(bin_num+1)
                                bin_halfwidth = 0.5*d_hist_yields[proc_type][proc_name].GetXaxis().GetBinWidth(bin_num+1)
                                
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
                    
                    syst_min_str, _ = pdgrnd.pdgRound(syst_min, 0)
                    syst_max_str, _ = pdgrnd.pdgRound(syst_max, 0)
                    syst_med_str, _ = pdgrnd.pdgRound(syst_med, 0)
                    
                    syst_min_str = str(syst_min_str) if (float(syst_min_str) >= 0.1) else "{<}0.1"
                    syst_max_str = str(syst_max_str) if (float(syst_max_str) >= 0.1) else "{<}0.1"
                    syst_med_str = str(syst_med_str) if (float(syst_med_str) >= 0.1) else "{<}0.1"
                    
                    # Set to the median if they differ by less than 1%
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
            ("toprule", "hline"),
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
            hist_obs = d_hist_yields["obs"]["obs"] if has_obs else None
            
            if has_obs :
                
                hist_ratio = d_hist_yields["obs"]["obs"].Clone("ratio")
            
            else :
                
                hist_ratio = d_hist_yields["bkg"]["total_bkg"].Clone("ratio")
            
            gr_bkg_unc = d_gr_uncs["bkg"]["total_bkg"]
            gr_ratio_unc = d_gr_uncs["bkg"]["total_bkg"].Clone("ratio_unc")
            
            for bin in range(hist_ratio.GetNbinsX()) :
                
                num = hist_ratio.GetBinContent(bin+1)
                num_unc = hist_ratio.GetBinError(bin+1)
                
                den = d_hist_yields["bkg"]["total_bkg"].GetBinContent(bin+1)
                den_unc_u = d_gr_uncs["bkg"]["total_bkg"].GetErrorYhigh(bin)
                den_unc_d = d_gr_uncs["bkg"]["total_bkg"].GetErrorYlow(bin)
                
                ratio = num/den if den else 0
                
                hist_ratio.SetBinContent(bin+1, ratio)
                hist_ratio.SetBinError(
                    bin+1,
                    ratio*num_unc/num if num else 0
                )
                
                gr_ratio_unc.SetPointY(bin, 1)
                gr_ratio_unc.SetPointEYhigh(bin, den_unc_u/den if den else 0)
                gr_ratio_unc.SetPointEYlow(bin, den_unc_d/den if den else 0)
            
            
            for hist in [hist_obs, hist_ratio] :
                
                if hist is None :
                    
                    continue
                
                hist.SetLineColor(1)
                hist.SetLineWidth(2)
                hist.SetMarkerColor(1)
                hist.SetMarkerSize(2)
                hist.SetMarkerStyle(20)
                hist.SetFillStyle(0)
                hist.SetOption("PE1")
            
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
            
            xmin = d_hist_yields["bkg"]["total_bkg"].GetXaxis().GetXmin()
            xmax = d_hist_yields["bkg"]["total_bkg"].GetXaxis().GetXmax()
            
            ymin = 1e-2
            ymax = 10**(round(numpy.log10(d_hist_yields["bkg"]["total_bkg"].GetMaximum()))+3)
            
            gr_unity = ROOT.TGraph()
            gr_unity.AddPoint(xmin, 1)
            gr_unity.AddPoint(xmax, 1)
            gr_unity.SetLineColor(1)
            gr_unity.SetLineWidth(1)
            gr_unity.GetHistogram().SetOption("L")
            
            outfilename = f"{outdir}/{infilename_noext}_plot.pdf"
            
            hist_dummy =  ROOT.TH1F(f"h1_dummy", "h1_dummy", nbins, 1, 1+nbins)
            
            l_hist_obs = [hist_obs]*int(has_obs)
            l_hist_bkg = [d_hist_yields["bkg"][_proc] for _proc in d_config["plot_bkgs"]]
            l_hist_sig = list(d_hist_yields["sig"].values())
            
            utils.root_plot1D_legacy(
                l_hist = l_hist_bkg,
                l_hist_overlay = l_hist_obs + l_hist_sig,
                ratio_num_den_pairs = [],
                l_graph_overlay = [gr_bkg_unc],
                l_ratio_hist_overlay = [hist_ratio] if has_obs else [hist_dummy],
                l_ratio_graph_overlay = [gr_unity, gr_ratio_unc],
                l_legend_order = l_hist_obs + [0]*int(has_obs) + l_hist_bkg + [gr_bkg_unc] + l_hist_sig,
                outfile = outfilename,
                xrange = (xmin, xmax),
                yrange = (ymin, ymax),
                logx = False,
                logy = d_config.get("logy", False),
                #xtitle = "Bin number",
                ytitle = "Events",
                #xtitle_ratio = "Bin number",
                xtitle_ratio = d_config["bin_header"],
                ytitle_ratio = "Data / Total bkg.",
                centertitlex = True, centertitley = True,
                centerlabelx = d_config.get("centerlabelx", False),
                centerlabely = False,
                gridx = False, gridy = False,
                ndivisionsx = (nbins, 1, 0),
                stackdrawopt = "",
                #ratiodrawopt = "E1",
                ndivisionsy_ratio = (2, 10, 0),
                legendpos = "UL",
                #legendtitle = "[m_{#tilde{#tau}}(250), m_{LSP}(1)]",
                legendncol = 2,
                #legendtextsize = 0.04,
                legendwidthscale = 1.8,
                legendheightscale = 5,
                no_xerror = True,
                CMSextraText = args.cmsextratext,
                lumiText = args.lumitext
            )
        
        
    return 0


if (__name__ == "__main__") :
    
    main()