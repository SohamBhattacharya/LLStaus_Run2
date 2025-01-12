#!/usr/bin/env -S python3 -u

import argparse
import numpy
import os
import sortedcontainers
import uproot

import ROOT
#import CombineHarvester.CombineTools.ch as ch

import utils.commonutils as cmut
import utils.utils as utils


def main() :
    
    # Argument parser
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        "--input",
        help = "Input root files",
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
    
    d_graph = sortedcontainers.SortedDict()
    
    l_wps = sortedcontainers.SortedList()
    l_dxys = sortedcontainers.SortedList()
    
    for ipoint, infilename in enumerate(args.input) :
        
        parsed_result = cmut.parse_string_regex(
            s = infilename,
            regexp = "ZMT_wp-(?P<wp>\w+)_dxy-gt-(?P<dxy>\w+)/",
        )
        
        wp = parsed_result["wp"]
        dxy = parsed_result["dxy"]
        
        wp_key = f"wp-{wp}"
        dxy_key = f"dxy-gt-{dxy}"
        
        wp = float(wp.replace("p", "."))
        dxy = float(dxy.replace("p", "."))
        
        if (wp_key not in l_wps) :
            
            l_wps.add(wp_key)
        
        if (dxy_key not in l_dxys) :
            
            l_dxys.add(dxy_key)
        
        with uproot.open(f"{infilename}:limit") as infile :
            
            a_deltaNLL = 2.0 * infile["deltaNLL"].array().to_numpy()
            a_sf = infile["SF"].array().to_numpy()
            
            npoint = len(a_sf)
            
            sf0 = a_sf[numpy.argmin(a_deltaNLL)]
            a_sf -= sf0
            
            sort_idx = numpy.argsort(a_sf)
            
            a_sf = a_sf[sort_idx]
            a_deltaNLL = a_deltaNLL[sort_idx]
            
            #print(sf0)
            #print(a_deltaNLL)
            #print(a_SF)
            
            gr = ROOT.TGraph(npoint, a_sf, a_deltaNLL)
            gr.SetName(f"{wp}_{dxy}")
            gr.SetTitle(f"{wp}_{dxy}")
            
            if (wp_key, dxy_key) not in d_graph :
                
                d_graph[(wp_key, dxy_key)] = gr
                d_graph[(wp_key, dxy_key)].SetName(f"gr_deltaNLL_{wp_key}_{dxy_key}")
                d_graph[(wp_key, dxy_key)].SetTitle(f"DisTau>{wp:0.2f}, d_{{xy}}>{dxy} cm")
            
            outfile.cd()
            gr.Write()
    
    
    outfile.Close()
    
    
    xmin = -0.5
    xmax = 0.5
    hist_dummy = ROOT.TH1F("hist_dummy", "", 10, xmin, xmax)
    
    for idxy, dxy_key in enumerate(l_dxys) :
        
        l_graph = []
        
        for iwp, wp_key in enumerate(l_wps) :
            
            gr = d_graph[(wp_key, dxy_key)]
            
            #for ipoint in range(0, gr.GetN()) :
            #    
            #    if (gr.GetPointY(ipoint) > 6) :
            #        
            #        print("Removing point {ipoint}")
            #        gr.RemovePoint(ipoint)
            
            color = utils.get_cms_colors(iwp)
            gr.SetLineColor(color)
            gr.SetMarkerColor(color)
            gr.SetMarkerSize(2)
            gr.SetMarkerStyle(20)
            gr.SetLineWidth(2)
            gr.SetLineStyle(7)
            gr.SetDrawOption("LP")
            gr.SetFillStyle(0)
            
            l_graph.append(gr)
        
        plotfile = f"{outdir}/deltaNLL_{dxy_key}.pdf"
        
        utils.root_plot1D_legacy(
            l_hist = [hist_dummy],
            l_graph_overlay = l_graph,
            gr_overlay_drawopt = "L",
            outfile = plotfile,
            xrange = (xmin, xmax),
            yrange = (-1, 10),
            #no_xerror = True,
            logx = False, logy = False,
            xtitle = "SF - #hat{SF}", ytitle = "-2#DeltaLL",
            centertitlex = True, centertitley = True,
            centerlabelx = False, centerlabely = False,
            gridx = True, gridy = True,
            ndivisionsx = (10, 5, 0),
            ndivisionsy = (11, 5, 0),
            stackdrawopt = "",
            legendpos = "UL",
            legendncol = 2,
            legendfillstyle = 1001,
            legendtextsize = 0.03,
            legendtitle = args.title,
            legendheightscale = 1.5, legendwidthscale = 1.9,
            CMSextraText = "Preliminary",
            lumiText = f"{args.era} (13 TeV)"
        )
        
        os.system(f"pdftoppm -png -r 600 -cropbox {plotfile} {plotfile}")
    
    
    return 0


if (__name__ == "__main__") :
    
    main()