#!/usr/bin/env python3

import argparse
#from traitlets import default
from ctypes import util
import functools
import os
import numpy
import operator
import sortedcontainers

import ROOT
#import CombineHarvester.CombineTools.ch as ch

import utils.commonutils as cmut
import utils.utils as utils

def main() :
    
    # Argument parser
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        "--jsons",
        help = "Input",
        type = str,
        nargs = "+",
        required = True,
    )
    
    parser.add_argument(
        "--output",
        help = "Output root file",
        type = str,
        required = True,
    )
    
    parser.add_argument(
        "--era",
        help = "Era string",
        type = str,
        required = True,
    )
    
    parser.add_argument(
        "--title",
        help = "Title",
        type = str,
        required = False,
        default = "",
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    outdir = os.path.dirname(args.output)
    
    if (len(outdir)) :
        
        os.system(f"mkdir -p {outdir}")
    
    outfile = ROOT.TFile.Open(args.output, "RECREATE")
    
    #cnv = ROOT.TCanvas("canvas", "canvas", 800, 800)
    
    gr2_sf = ROOT.TGraph2DAsymmErrors()
    gr2_sf.SetName("gr_SFs")
    gr2_sf.SetTitle("DisTau Scale Factors;WP;d_{xy} cut;SF")
    
    #gr2_sf.GetXaxis().SetTitle("WP")
    #gr2_sf.GetYaxis().SetTitle("d_{xy} cut")
    #gr2_sf.GetZaxis().SetTitle("SF")
    
    h2 = ROOT.TH2F("h2_SFs", "DisTau Scale Factors", 500, 0.6, 1.1, 100, 0.0, 0.1)
    
    h2.GetXaxis().SetTitle("WP")
    h2.GetYaxis().SetTitle("d_{xy} cut")
    h2.GetZaxis().SetTitle("SF")
    
    d_gr_wps = sortedcontainers.SortedDict()
    d_gr_dxys = sortedcontainers.SortedDict()
    
    for ipoint, injson in enumerate(args.jsons) :
        
        fit_results = cmut.load_config(injson)
        
        if "120.0" not in fit_results :
            continue
        
        fit_results = fit_results["120.0"]
        print(fit_results)
        
        parsed_result = cmut.parse_string_regex(
            s = injson,
            regexp = "fit_result_ZMT_wp-(?P<wp>\w+)_dxy-gt-(?P<dxy>\w+).json",
        )
        
        wp = parsed_result["wp"]
        dxy = parsed_result["dxy"]
        
        wp_key = f"wp-{wp}"
        dxy_key = f"dxy-gt-{dxy}"
        
        wp = float(wp.replace("p", "."))
        dxy = float(dxy.replace("p", "."))
        
        if wp_key not in d_gr_wps :
            
            d_gr_wps[wp_key] = ROOT.TGraphAsymmErrors()
            d_gr_wps[wp_key].SetName(f"gr_SFs_{wp_key}")
            d_gr_wps[wp_key].SetTitle(f"DisTau > {wp:0.2f}")
        
        if dxy_key not in d_gr_dxys :
            
            d_gr_dxys[dxy_key] = ROOT.TGraphAsymmErrors()
            d_gr_dxys[dxy_key].SetName(f"gr_SFs_{dxy_key}")
            d_gr_dxys[dxy_key].SetTitle(f"Jet d_{{xy}} > {dxy:0.2f} cm")
        
        #wp = float(wp.replace("p", "."))
        #dxy = float(dxy.replace("p", "."))
        
        sf = fit_results["exp0"]
        sf_upr = fit_results["exp+1"]
        sf_lwr = fit_results["exp-1"]
        
        sf_err_lwr = abs(sf - sf_lwr)
        sf_err_upr = abs(sf - sf_upr)
        
        print(wp, dxy, sf)
        print(injson, sf, sf_upr, sf_lwr)
        
        gr2_sf.SetPoint(ipoint, wp, dxy, sf)
        gr2_sf.SetPointError(ipoint, 0.0, 0.0, 0.0, 0.0, sf_err_lwr, sf_err_upr)
        
        ipoint_tmp = d_gr_wps[wp_key].GetN()
        print(wp_key, ipoint_tmp)
        d_gr_wps[wp_key].SetPoint(ipoint_tmp, dxy, sf)
        d_gr_wps[wp_key].SetPointEYhigh(ipoint_tmp, sf_err_upr)
        d_gr_wps[wp_key].SetPointEYlow(ipoint_tmp, sf_err_lwr)
        
        ipoint_tmp = d_gr_dxys[dxy_key].GetN()
        print(dxy_key, ipoint_tmp)
        d_gr_dxys[dxy_key].SetPoint(ipoint_tmp, wp, sf)
        d_gr_dxys[dxy_key].SetPointEYhigh(ipoint_tmp, sf_err_upr)
        d_gr_dxys[dxy_key].SetPointEYlow(ipoint_tmp, sf_err_lwr)
        
        h2.Fill(wp, dxy, sf)
        h2.SetBinError(
            h2.GetXaxis().FindBin(wp),
            h2.GetYaxis().FindBin(dxy),
            0.5*(sf_err_lwr+sf_err_upr)
        )
    
    #h2_fromgr = gr2_sf.GetHistogram("empty").Clone()
    #h2_fromgr2_sf.SetName("SFs_fromgr")
    
    gr2_sf.SetNpx(20)
    gr2_sf.SetNpy(20)
    
    h2_interp = gr2_sf.GetHistogram().Clone()
    h2_interp.SetName("SFs_interpolated")
    
    outfile.cd()
    
    gr2_sf.Write()
    h2.Write()
    h2_interp.Write()
    #h2_fromgr2_sf.Write()
    
    for igr, (key, gr) in enumerate(d_gr_wps.items()) :
        
        for ipoint in range(0, gr.GetN()) :
            
            gr.SetPointX(
                ipoint,
                gr.GetPointX(ipoint)+(0.01/len(d_gr_wps)*igr)
            )
        
        color = utils.get_cms_colors(igr)
        gr.SetLineColor(color)
        gr.SetMarkerColor(color)
        gr.SetMarkerSize(2)
        gr.SetMarkerStyle(20)
        gr.SetLineWidth(2)
        gr.SetDrawOption("PE1")
        gr.SetFillStyle(0)
        #gr.SetFillColor(0)
        
        gr.Write()
    
    for igr, (key, gr) in enumerate(d_gr_dxys.items()) :
        
        for ipoint in range(0, gr.GetN()) :
            
            gr.SetPointX(
                ipoint,
                gr.GetPointX(ipoint)+(0.01/len(d_gr_wps)*igr)
            )
        
        color = utils.get_cms_colors(igr)
        gr.SetLineColor(color)
        gr.SetMarkerColor(color)
        gr.SetMarkerSize(2)
        gr.SetMarkerStyle(20)
        gr.SetLineWidth(2)
        gr.SetDrawOption("PE1")
        gr.SetFillStyle(0)
        
        gr.Write()
    
    outfile.Close()
    
    min_dxy = 0.0
    max_dxy = 0.1
    hist_vs_dxy = ROOT.TH1F("hist_vs_dxy", "", 5, min_dxy, max_dxy)
    plotfile = f"{outdir}/DisTauSF_vs_dxy.pdf"
    
    utils.root_plot1D_legacy(
        l_hist = [hist_vs_dxy],
        l_graph_overlay = list(d_gr_wps.values()),
        outfile = plotfile,
        xrange = (min_dxy, max_dxy),
        yrange = (0, 2.5),
        #no_xerror = True,
        logx = False, logy = False,
        xtitle = "Jet d_{xy} threshold [cm] group", ytitle = "SF",
        centertitlex = True, centertitley = True,
        centerlabelx = False, centerlabely = False,
        gridx = True, gridy = True,
        ndivisionsx = (10, 5, 0),
        #ndivisionsy = (4, 5, 0),
        stackdrawopt = "",
        legendpos = "UL",
        legendncol = 2,
        legendtextsize = 0.03,
        legendtitle = args.title,
        legendheightscale = 1.5, legendwidthscale = 1.5,
        #CMSextraText = "Preliminary",
        CMSextraText = "Private Work",
        lumiText = f"{args.era} (13 TeV)"
    )
    
    #os.system(f"pdftoppm -png -r 600 -cropbox {plotfile} {plotfile}")
    
    min_wp = 0.6
    max_wp = 1.1
    hist_vs_wp = ROOT.TH1F("hist_vs_wp", "", 5, min_wp, max_wp)
    plotfile  = f"{outdir}/DisTauSF_vs_wp.pdf"
    
    utils.root_plot1D_legacy(
        l_hist = [hist_vs_wp],
        l_graph_overlay = list(d_gr_dxys.values()),
        outfile = f"{outdir}/DisTauSF_vs_wp.pdf",
        xrange = (min_wp, max_wp),
        yrange = (0, 2.5),
        #no_xerror = True,
        logx = False, logy = False,
        xtitle = "DisTau threshold group", ytitle = "SF",
        centertitlex = True, centertitley = True,
        centerlabelx = False, centerlabely = False,
        gridx = True, gridy = True,
        ndivisionsx = (10, 5, 0),
        #ndivisionsy = (4, 5, 0),
        stackdrawopt = "",
        legendpos = "UL",
        legendncol = 2,
        legendtextsize = 0.03,
        legendtitle = args.title,
        legendheightscale = 1.5, legendwidthscale = 1.5,
        #CMSextraText = "Preliminary",
        CMSextraText = "Private Work",
        lumiText = f"{args.era} (13 TeV)"
    )
    
    #os.system(f"pdftoppm -png -r 600 -cropbox {plotfile} {plotfile}")
    
    return 0


if (__name__ == "__main__") :
    
    main()