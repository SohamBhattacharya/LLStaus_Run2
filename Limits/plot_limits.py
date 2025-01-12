#!/usr/bin/env python3

import argparse
import pathlib
import numpy
import os
import sortedcontainers

import ROOT
ROOT.gROOT.SetBatch(True)

import utils.commonutils as cmut


def getContourListFromHist(hist, contourVal, nSmooth = 0) :
    
    #print "*****", hist.GetNbinsX(), hist.GetNbinsY()
    
    #current_pad = ROOT.gROOT.GetSelectedPad()
    
    if (nSmooth > 0) :

        #hist.Smooth(1, "k5a")
        #hist.Smooth(nSmooth, "k5b")

        # Works better
        for iSmooth in range(0, nSmooth) :

            #hist.Smooth(1, "k5b")
            hist.Smooth(1, "kba")
            hist.Smooth(1, "k5a")
    
    canvas_temp = ROOT.TCanvas("canvas_temp", "canvas_temp")
    canvas_temp.cd()
    
    l_contourVal = [contourVal]
    
    hist_temp = hist.Clone()
    ROOT.SetOwnership(hist_temp, 0)
    
    # Set contours
    hist_temp.SetContour(len(l_contourVal), numpy.array(l_contourVal, dtype = numpy.float64))
    hist_temp.Draw("CONT Z LIST")
    canvas_temp.Update()
    
    # Get contours
    ROOT.gROOT.GetListOfSpecials().Print()
    l_contour = ROOT.gROOT.GetListOfSpecials().FindObject("contours").Clone()
    l_contour = l_contour.At(0)
    
    #if (not l_contour.GetEntries()) :
    #    
    #    print("Error: No contours found.")
    #    exit(1)
    
    canvas_temp.Clear()
    canvas_temp.Close()
    
    # Switch back to the original pad
    #current_pad.cd()
    
    return l_contour


def getContourListFromGraph(graph, contourVal, nSmooth = 0) :
    
    return getContourListFromHist(graph.GetHistogram().Clone(), contourVal, nSmooth)


def get_exclusion(
    jsons,
    json_key,
    name,
    axis_label,
) :
    
    gr = ROOT.TGraph2D()
    gr.SetName(f"gr_{name}")
    
    
    h2 = ROOT.TH2F(f"h2_{name}", f"h2_{name}", 60, 0, 600, int(2e5), 0.01, 2000)
    h2.GetXaxis().SetTitle("m(#tilde{#tau}) [GeV]")
    h2.GetYaxis().SetTitle("c#tau_{0} [mm]")
    h2.GetZaxis().SetTitle(axis_label)
    
    for ipoint, injson in enumerate(jsons) :
        
        limits = cmut.load_config(injson)["120.0"]
        print(limits)
        
        sig_info = cmut.parse_stau_samplestring(
            s = injson,
            regexp = "SMS-TStauStau_MStau-(?P<mstau>\d+)_ctau-(?P<ctau>\w+)_mLSP-(?P<mlsp>\d+)",
        )
        
        mstau = sig_info["mstau"]
        ctau = sig_info["ctau"]
        ctau = float(ctau[0: -2].replace("p", "."))
        
        ul = limits[json_key] if json_key in limits else 0
        ul = numpy.clip(ul, 0.0, 2.0)
        print(injson, json_key, mstau, ctau, ul)
        
        gr.SetPoint(ipoint, mstau, ctau, ul)
        h2.Fill(mstau, ctau, ul)
    
    gr.SetNpx(500)
    gr.SetNpy(500)
    
    h2_interp = gr.GetHistogram().Clone()
    h2_interp.SetName(f"{h2.GetName()}_interp")
    h2_interp.GetXaxis().SetTitle(h2.GetXaxis().GetTitle())
    h2_interp.GetYaxis().SetTitle(h2.GetYaxis().GetTitle())
    h2_interp.GetZaxis().SetTitle(h2.GetZaxis().GetTitle())
    
    contours = getContourListFromHist(
        hist = h2_interp.Clone(),
        contourVal = 1.0,
        nSmooth = 0,
    )
    
    result = sortedcontainers.SortedDict()
    
    result["gr"] = gr
    result["h2"] = h2
    result["h2_interp"] = h2_interp
    
    if (len(contours)) :
        
        cnt = contours[0]
        cnt.SetName(f"contour_{name}")
        result["contour"] = cnt
    
    return result


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
    
    # Parse arguments
    args = parser.parse_args()
    
    outdir = os.path.dirname(args.output).strip()
    
    if (len(outdir)) :
        
        os.system(f"mkdir -p {outdir}")
    
    outfile = ROOT.TFile.Open(args.output, "RECREATE")
    
    results = sortedcontainers.SortedDict()
    
    results["obs"] = get_exclusion(
        jsons = args.jsons,
        json_key = "obs",
        name = "obs",
        axis_label = "Observed upper limit",
    )
    
    results["exp"] = get_exclusion(
        jsons = args.jsons,
        json_key = "exp0",
        name = "exp",
        axis_label = "Expected upper limit",
    )
    
    results["exp_p1"] = get_exclusion(
        jsons = args.jsons,
        json_key = "exp+1",
        name = "exp_p1",
        axis_label = "Expected upper limit (+1#sigma)",
    )
    
    results["exp_m1"] = get_exclusion(
        jsons = args.jsons,
        json_key = "exp-1",
        name = "exp_m1",
        axis_label = "Expected upper limit (-1#sigma)",
    )
    
    outfile.cd()
    
    print(f"Writing to: {args.output}")
    
    for result_key, result in results.items() :
        
        outfile.cd()
        outfile.mkdir(result_key)
        outfile.cd(result_key)
        
        for key, val in result.items() :
            
            print(f"Writing [{result_key}] [{key}]")
            val.Write()
    
    return 0


if (__name__ == "__main__") :
    
    main()