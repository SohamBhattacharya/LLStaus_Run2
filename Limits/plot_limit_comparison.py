#!/usr/bin/env python3

import argparse
import ctypes
import cmsstyle as CMS
import numpy
import os
import ROOT

from matplotlib import contour
from sympy import root
import utils.commonutils as cmut
import utils.utils as utils


def main() :
    
    # Argument parser
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        "--outdir",
        help = "Output directory",
        type = str,
        required = True,
    )
    
    parser.add_argument(
        "--extratext",
        help = "Extra text to be printed on canvas",
        type = str,
        required = False,
    )
    
    parser.add_argument(
        "--cmsextratext",
        help = "Will set the CMS extra text; can be empty string \"\"",
        type = str,
        required = False,
        default = "Preliminary",
    )
    
    args = parser.parse_args()
    
    os.system(f"mkdir -p {args.outdir}")
    
    d_contour_info = {}
    
    d_contour_info["all_bins"] = {
        "fname": "results/limits_signal_v9_w-distau-sfs/llstau_maximally-mixed/channels_BRT2/eras_all/limits/limits.root",
        "contour": "exp/contour_exp",
        "h2_ul": "exp/h2_ul_exp_interp",
        "label": "All bins",
        #"color": utils.get_cms_colors(len(d_contour_info))
        "color": ROOT.kBlack,
        "linestyle": ROOT.kSolid
    }
    d_contour_info["wo_bin1"] = {
        "fname": "results/limits_signal_v9_w-distau-sfs/llstau_maximally-mixed/channels_BRT2/eras_all/limits_disable_bin3_BRT2/limits_disable_bin3_BRT2.root",
        "contour": "exp/contour_exp",
        "h2_ul": "exp/h2_ul_exp_interp",
        "label": "Without bin 1",
        #"color": utils.get_cms_colors(len(d_contour_info))-1
        "color": ROOT.kBlack,
        "linestyle": ROOT.kDashed
    }
    d_contour_info["wo_bin2"] = {
        "fname": "results/limits_signal_v9_w-distau-sfs/llstau_maximally-mixed/channels_BRT2/eras_all/limits_disable_bin4_BRT2/limits_disable_bin4_BRT2.root",
        "contour": "exp/contour_exp",
        "h2_ul": "exp/h2_ul_exp_interp",
        "label": "Without bin 2",
        #"color": utils.get_cms_colors(len(d_contour_info))-1
        "color": ROOT.kBlack,
        "linestyle": ROOT.kDashed
    }
    d_contour_info["wo_bin3"] = {
        "fname": "results/limits_signal_v9_w-distau-sfs/llstau_maximally-mixed/channels_BRT2/eras_all/limits_disable_bin5_BRT2/limits_disable_bin5_BRT2.root",
        "contour": "exp/contour_exp",
        "h2_ul": "exp/h2_ul_exp_interp",
        "label": "Without bin 3",
        #"color": utils.get_cms_colors(len(d_contour_info))-1
        "color": ROOT.kBlack,
        "linestyle": ROOT.kDashed
    }
    d_contour_info["wo_bin4"] = {
        "fname": "results/limits_signal_v9_w-distau-sfs/llstau_maximally-mixed/channels_BRT2/eras_all/limits_disable_bin6_BRT2/limits_disable_bin6_BRT2.root",
        "contour": "exp/contour_exp",
        "h2_ul": "exp/h2_ul_exp_interp",
        "label": "Without bin 4",
        #"color": utils.get_cms_colors(len(d_contour_info))-1
        "color": ROOT.kBlack,
        "linestyle": ROOT.kDashed
    }
    d_contour_info["wo_bin5"] = {
        "fname": "results/limits_signal_v9_w-distau-sfs/llstau_maximally-mixed/channels_BRT2/eras_all/limits_disable_bin7_BRT2/limits_disable_bin7_BRT2.root",
        "contour": "exp/contour_exp",
        "h2_ul": "exp/h2_ul_exp_interp",
        "label": "Without bin 5",
        #"color": utils.get_cms_colors(len(d_contour_info))-1
        "color": ROOT.kBlack,
        "linestyle": ROOT.kDashed
    }
    d_contour_info["wo_bin6"] = {
        "fname": "results/limits_signal_v9_w-distau-sfs/llstau_maximally-mixed/channels_BRT2/eras_all/limits_disable_bin8_BRT2/limits_disable_bin8_BRT2.root",
        "contour": "exp/contour_exp",
        "h2_ul": "exp/h2_ul_exp_interp",
        "label": "Without bin 6",
        #"color": utils.get_cms_colors(len(d_contour_info))-1
        "color": ROOT.kBlack,
        "linestyle": ROOT.kDashed
    }
    d_contour_info["wo_bin7"] = {
        "fname": "results/limits_signal_v9_w-distau-sfs/llstau_maximally-mixed/channels_BRT2/eras_all/limits_disable_bin9_BRT2/limits_disable_bin9_BRT2.root",
        "contour": "exp/contour_exp",
        "h2_ul": "exp/h2_ul_exp_interp",
        "label": "Without bin 7",
        #"color": utils.get_cms_colors(len(d_contour_info))-1
        "color": ROOT.kBlack,
        "linestyle": ROOT.kDashed
    }
    d_contour_info["wo_bin8"] = {
        "fname": "results/limits_signal_v9_w-distau-sfs/llstau_maximally-mixed/channels_BRT2/eras_all/limits_disable_bin10_BRT2/limits_disable_bin10_BRT2.root",
        "contour": "exp/contour_exp",
        "h2_ul": "exp/h2_ul_exp_interp",
        "label": "Without bin 8",
        #"color": utils.get_cms_colors(len(d_contour_info))-1
        "color": ROOT.kBlack,
        "linestyle": ROOT.kDashed
    }
    
    
    xmin = 90
    xmax = 600
    ymin = 1
    ymax = 10**5
    #zmin = -0.1
    #zmax = 0.5
    #zdiv = 0.1
    zmin = -10
    zmax = 50
    zdiv = 10
    zaxis_ncontours = round((zmax - zmin) / zdiv) #* 2
    
    #CMS.setCMSStyle()
    #CMS.SetExtraText(args.cmsextratext)
    #CMS.SetLumi("")
    #CMS.SetEnergy("")
    #CMS.SetCMSPalette()
    #CMS.ResetAdditionalInfo()
    #
    #canvas_name = "limits-vs-ctau0-mstau"
    #canvas = CMS.cmsCanvas(
    #    canvName = canvas_name,
    #    x_min = xmin, x_max = xmax,
    #    y_min = ymin, y_max = ymax,
    #    nameXaxis = "m_{#tilde{#tau}} [GeV]",
    #    #nameYaxis = "c#tau_{0}(#tilde{#tau}) [mm]",
    #    nameYaxis = "c#tau_{0} [mm]",
    #    #square = CMS.kSquare,
    #    square = False,
    #    iPos = 0,
    #    extraSpace = 0.03,
    #    with_z_axis = True,
    #    scaleLumi = None,
    #    yTitOffset = 1.05
    #)
    #
    #canvas.SetRightMargin(0.175)
    #CMS.SetLumi(137.62, round_lumi = True)
    #CMS.SetEnergy(13)
    #CMS.CMS_lumi(canvas, iPosX = 0)
    #CMS.UpdatePad(canvas)
    #
    #h2_dummy = ROOT.TH2F("h2_dummy", "h2_dummy", 1, xmin, xmax, 1, ymin, ymax)
    #
    #zmin = 10**-1
    #zmax = 10**3
    #
    #h2_dummy.GetZaxis().SetRangeUser(zmin, zmax)
    #h2_dummy.GetZaxis().SetTitle("95% CL upper limit on cross section [fb]")
    #h2_dummy.GetZaxis().SetTitleSize(0.05)
    #h2_dummy.GetZaxis().SetTitleOffset(1.2)
    #h2_dummy.GetZaxis().SetLabelSize(0.05)
    #
    #legend = CMS.cmsLeg(canvas.GetLeftMargin(), 0.65, 1-canvas.GetRightMargin(), 1-canvas.GetTopMargin(), textSize = 0.0425, columns = 2)
    #legend.SetFillColor(0)
    #legend.SetFillStyle(1001)
    #legend.SetBorderSize(1)
    #legend.SetTextAlign(22)
    
    l_draw_objects = []
    l_bins_contour_ref = []
    
    contour_key_ref = "all_bins"
    h2_ul_ref = None
    contour_ref = None
    
    for icnt, (contour_key, contour_info) in enumerate(d_contour_info.items()) :
        
        label = contour_info["label"]
        color = contour_info["color"]
        linestyle = contour_info["linestyle"]
        rootfile = ROOT.TFile.Open(contour_info["fname"], "READ")
        
        contour_info_ref = d_contour_info[contour_key_ref]
        
        contour = rootfile.Get(contour_info["contour"])
        #contour.SetDirectory(0)
        
        h2_ul = rootfile.Get(contour_info["h2_ul"])
        h2_ul.SetDirectory(0)
        
        if not h2_ul_ref:
            
            h2_ul_ref = h2_ul.Clone()
            h2_ul_ref.SetDirectory(0)
            contour_ref = contour.Clone()
            continue
        
        h2_ul.Add(h2_ul_ref, -1)
        h2_ul.Divide(h2_ul_ref)
        h2_ul.Scale(100)
        
        #palette = ROOT.Experimental.TPalette(ROOT.Experimental.TPalette.kDiscrete, [[-0.1, ROOT.kYellow], [0.5, ROOT.kGreen]])
        
        CMS.setCMSStyle()
        CMS.SetExtraText(args.cmsextratext)
        CMS.SetLumi("")
        CMS.SetEnergy("")
        CMS.SetCMSPalette()
        CMS.ResetAdditionalInfo()
        #CMS.getCMSStyle().SetPalette(ROOT.kRainBow)
        CMS.getCMSStyle().SetPalette(ROOT.kLightTemperature)
        #CMS.getCMSStyle().SetPalette(palette)
        CMS.getCMSStyle().SetNumberContours(zaxis_ncontours)
        
        canvas_name = f"limits-vs-ctau0-mstau_compare_{contour_key}.{contour_key_ref}"
        canvas = CMS.cmsCanvas(
            canvName = canvas_name,
            x_min = xmin, x_max = xmax,
            y_min = ymin, y_max = ymax,
            nameXaxis = "m_{#tilde{#tau}} [GeV]",
            nameYaxis = "c#tau_{0} [mm]",
            square = False,
            iPos = 0,
            extraSpace = 0.03,
            with_z_axis = True,
            scaleLumi = None,
            yTitOffset = 1.05
        )
        
        canvas.SetRightMargin(0.175)
        CMS.SetLumi(137.62, round_lumi = True)
        CMS.SetEnergy(13)
        CMS.CMS_lumi(canvas, iPosX = 0)
        CMS.UpdatePad(canvas)
        
        h2_ul.GetZaxis().SetRangeUser(zmin, zmax)
        h2_ul.GetZaxis().CenterTitle(True)
        h2_ul.GetZaxis().SetTitle(f"(UL_{{{label}}} #minus UL_{{{contour_info_ref['label']}}}) / UL_{{{contour_info_ref['label']}}} [%]")
        h2_ul.GetZaxis().SetTitleSize(0.05)
        h2_ul.GetZaxis().SetTitleOffset(1.2)
        h2_ul.GetZaxis().SetLabelSize(0.05)
        h2_ul.GetZaxis().SetLabelOffset(0.01)
        
        legend = CMS.cmsLeg(canvas.GetLeftMargin(), 0.6, 1-canvas.GetRightMargin(), 1-canvas.GetTopMargin(), textSize = 0.0425, columns = 2)
        legend.SetFillColor(0)
        legend.SetFillStyle(1001)
        legend.SetBorderSize(1)
        legend.SetTextAlign(22)
        legend.SetMargin(0.4)
        #legend.SetColumnSeparation(0.5)
        
        CMS.cmsDraw(h2_ul, "colz")
        CMS.cmsDraw(contour_ref, style = "sameL", lstyle = contour_info_ref["linestyle"], lcolor = contour_info_ref["color"], lwidth = 3)
        CMS.cmsDraw(contour, style = "sameL", lstyle = linestyle, lcolor = color, lwidth = 3)
        
        legend.AddEntry(contour_ref, f"#splitline{{Expected}}{{({contour_info_ref['label']})}}", "L")
        legend.AddEntry(contour, f"#splitline{{Expected}}{{({label})}}", "L")
        
        label = (
            "pp#rightarrow#tilde{#tau}#bar{#tilde{#tau}} , "
            "#tilde{#tau}#rightarrow#tau#tilde{G} , "
            "m_{#tilde{G}} = 1 GeV"
        )
        
        if args.extratext :
            
            label = f"#splitline{{{args.extratext}}}{{{label}}}"
        
        label = f"#lower[0.2]{{{label}}}"
        
        CMS.cmsHeader(
            leg = legend,
            legTitle = label,
            textAlign = 22,
            textSize = 0.0425,
            textFont = 62,
            #textColor = ROOT.kBlack,
            isToRemove = True
        )
        
        CMS.GetcmsCanvasHist(canvas).GetXaxis().SetNdivisions(6, 5, 0)
        
        canvas.SetLogy()
        legend.Draw()
        canvas.Update()
        
        #canvas.SetLogz()
        canvas_outfile = f"{args.outdir}/{canvas_name}"
        CMS.SaveCanvas(canvas, f"{canvas_outfile}.pdf")
        cmut.pdf_to_png(f"{canvas_outfile}.pdf")
        
        #nbinsx = h2_ul.GetNbinsX()
        #nbinsy = h2_ul.GetNbinsY()
        #
        #firstbinx = h2_ul.GetXaxis().FindBin(90)
        #firstbiny = h2_ul.GetXaxis().FindBin(1)
        #
        #binx = ctypes.c_int(-1)
        #biny = ctypes.c_int(-1)
        #
        #if not icnt :
        #    
        #    l_bins_contour_ref = [(_x, _y) for _x in range(1, nbinsx+1) for _y in range(1, nbinsy+1) if h2_ul.GetBinContent(_x, _y) < 1]
        #
        ##l_uls = [h2_ul.GetBinContent(_x, _y) for _x in range(1, nbinsx+1) for _y in range(1, nbinsy+1) if h2_ul.GetBinContent(_x, _y) < 1]
        #l_uls = [h2_ul.GetBinContent(*_bin) for _bin in l_bins_contour_ref]
        #
        ## Min
        ##minz = min(l_uls)
        #minz = h2_ul.GetMinimum()
        #h2_ul.GetBinWithContent2(minz, binx, biny, firstbinx, nbinsx+1, firstbiny, nbinsy+1, 0.01)
        #
        #minx = h2_ul.GetXaxis().GetBinCenter(binx.value)
        #miny = h2_ul.GetYaxis().GetBinCenter(biny.value)
        #
        #g1_min = ROOT.TGraph()
        #g1_min.SetName(f"g1_min_{icnt}")
        #g1_min.AddPoint(minx, miny)
        #l_draw_objects.append(g1_min)
        #
        ## Max
        #maxz = max(l_uls)
        #h2_ul.GetBinWithContent2(maxz, binx, biny, firstbinx, nbinsx+1, firstbiny, nbinsy+1, 0.01)
        #
        #maxx = h2_ul.GetXaxis().GetBinCenter(binx.value)
        #maxy = h2_ul.GetYaxis().GetBinCenter(biny.value)
        #
        #g1_max = ROOT.TGraph()
        #g1_max.SetName(f"g1_max_{icnt}")
        #g1_max.AddPoint(maxx, maxy)
        #l_draw_objects.append(g1_max)
        #
        ## Mean
        #meanz = numpy.mean(l_uls)
        #h2_ul.GetBinWithContent2(meanz, binx, biny, firstbinx, nbinsx+1, firstbiny, nbinsy+1, 0.01)
        #
        #meanx = h2_ul.GetXaxis().GetBinCenter(binx.value)
        #meany = h2_ul.GetYaxis().GetBinCenter(biny.value)
        #
        #g1_mean = ROOT.TGraph()
        #g1_mean.SetName(f"g1_mean_{icnt}")
        #g1_mean.AddPoint(meanx, meany)
        #l_draw_objects.append(g1_mean)
        #
        #ln_min_mean = ROOT.TLine(minx, miny, meanx, meany)
        #l_draw_objects.append(ln_min_mean)
        #
        #ln_mean_max = ROOT.TLine(meanx, meany, maxx, maxy)
        #l_draw_objects.append(ln_mean_max)
        #
        #print(label, h2_ul.GetMinimum(), h2_ul.GetMean(3), binx, biny, minx, miny)
        #
        #contour.SetLineColor(color)
        #contour.SetLineWidth(3)
        #contour.SetLineStyle(ROOT.kDashed)
        ##CMS.cmsDraw(contour, style = "sameL", lstyle = ROOT.kDashed, lcolor = color, lwidth = 3)
        #CMS.cmsDraw(g1_min, style = "sameP", mcolor = color, marker = 53, msize = 1)
        #CMS.cmsDraw(g1_mean, style = "sameP", mcolor = color, marker = 54, msize = 1)
        #CMS.cmsDraw(g1_max, style = "sameP", mcolor = color, marker = 55, msize = 1)
        #CMS.cmsDrawLine(ln_min_mean, lstyle = ROOT.kDotted, lcolor = color, lwidth = 1)
        #CMS.cmsDrawLine(ln_mean_max, lstyle = ROOT.kDotted, lcolor = color, lwidth = 1)
        #
        #legend.AddEntry(contour, label, "L")
        
        rootfile.Close()
    
    #CMS.GetcmsCanvasHist(canvas).GetXaxis().SetNdivisions(6, 5, 0)
    #
    #canvas.SetLogy()
    #
    ##label = (
    ##    "pp#rightarrow#tilde{#tau}#bar{#tilde{#tau}} , "
    ##    "#tilde{#tau}#rightarrow#tau#tilde{G} , "
    ##    "m_{#tilde{G}} = 1 GeV"
    ##)
    ##
    ##if args.extratext :
    ##    
    ##    label = f"#splitline{{{args.extratext}}}{{{label}}}"
    ##
    ##label = f"#lower[0.2]{{{label}}}"
    ##
    ##CMS.cmsHeader(
    ##    leg = legend,
    ##    legTitle = label,
    ##    textAlign = 22,
    ##    textSize = 0.0425,
    ##    textFont = 62,
    ##    #textColor = ROOT.kBlack,
    ##    isToRemove = True
    ##)
    #
    ##canvas.RedrawAxis()
    #legend.Draw()
    #canvas.Update()
    #
    #canvas.SetLogz()
    #
    ##outfile.cd()
    ##canvas.Write()
    #
    #canvas_outfile = f"{args.outdir}/{canvas_name}"
    #CMS.SaveCanvas(canvas, f"{canvas_outfile}.pdf")
    #cmut.pdf_to_png(f"{canvas_outfile}.pdf")
    
    return 0

if __name__ == "__main__" :
    
    main()