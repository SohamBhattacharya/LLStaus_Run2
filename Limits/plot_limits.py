#!/usr/bin/env python3

import argparse
import cmsstyle as CMS
import numpy
import os
import sortedcontainers

import ROOT
ROOT.gROOT.SetBatch(True)

import utils.commonutils as cmut


def get_contours_from_hist(hist, contour_value, nsmooth = 0) :
    
    if (nsmooth > 0) :

        #hist.Smooth(1, "k5a")
        #hist.Smooth(nsmooth, "k5b")

        # Works better
        for ismooth in range(0, nsmooth) :

            #hist.Smooth(1, "k5b")
            hist.Smooth(1, "kba")
            hist.Smooth(1, "k5a")
    
    canvas_temp = ROOT.TCanvas("canvas_temp", "canvas_temp")
    canvas_temp.cd()
    
    l_contours_values = [contour_value]
    
    hist_tmp = hist.Clone()
    ROOT.SetOwnership(hist_tmp, 0)
    
    # Set contours
    hist_tmp.SetContour(len(l_contours_values), numpy.array(l_contours_values, dtype = numpy.float64))
    hist_tmp.Draw("CONT Z LIST")
    canvas_temp.Update()
    
    # Get contours
    ROOT.gROOT.GetListOfSpecials().Print()
    l_contours = ROOT.gROOT.GetListOfSpecials().FindObject("contours").Clone()
    l_contours = l_contours.At(0)
    
    #if (not l_contours.GetEntries()) :
    #    
    #    print("Error: No contours found.")
    #    exit(1)
    
    canvas_temp.Clear()
    canvas_temp.Close()
    
    return l_contours


def get_contours_from_graph(graph, contour_value, nsmooth = 0) :
    
    return get_contours_from_hist(graph.GetHistogram().Clone(), contour_value, nsmooth)


def get_exclusion(
    d_limits,
    d_xsecs,
    limit_key,
    name,
    axis_label,
) :
    g2_ul = ROOT.TGraph2D()
    g2_ul.SetName(f"g2_ul_{name}")
    
    h2_ul = ROOT.TH2F(f"h2_ul_{name}", f"h2_ul_{name}", 100, 0, 1000, int(2e5), 0.001, 2000)
    h2_ul.GetXaxis().SetTitle("m(#tilde{#tau}) [GeV]")
    h2_ul.GetYaxis().SetTitle("c#tau_{0} [mm]")
    h2_ul.GetZaxis().SetTitle(axis_label)
    
    g2_xsecul = ROOT.TGraph2D()
    g2_xsecul.SetName(f"g2_xsecul_{name}")
    
    h2_xsecul = ROOT.TH2F(f"h2_xsecul_{name}", f"h2_xsecul_{name}", 100, 0, 1000, int(2e5), 0.001, 2000)
    h2_xsecul.GetXaxis().SetTitle("m(#tilde{#tau}) [GeV]")
    h2_xsecul.GetYaxis().SetTitle("c#tau_{0} [mm]")
    h2_xsecul.GetZaxis().SetTitle(f"{axis_label} on cross section [fb]")
    
    for ipoint, ((mstau, ctau), limits) in enumerate(d_limits.items()) :
        
        if limit_key not in limits :
            
            print(f"Warning: could not find limit {limit_key} for mstau={mstau}, ctau={ctau}. Skipping")
            continue
        
        ul = limits[limit_key] # if limit_key in limits else 0
        xsec = d_xsecs[mstau]["nom"]
        xsecul = xsec*ul
        print(limit_key, mstau, ctau, ul, xsecul)
        
        g2_ul.SetPoint(ipoint, mstau, ctau, ul)
        h2_ul.Fill(mstau, ctau, ul)
        
        g2_xsecul.SetPoint(ipoint, mstau, ctau, xsecul)
        h2_xsecul.Fill(mstau, ctau, xsecul)
    
    #g2_ul.GetHistogram().GetXaxis().SetRangeUser(0, 1000)
    #g2_xsecul.GetHistogram().GetXaxis().SetRangeUser(0, 1000)
    
    g2_ul.SetNpx(500)
    g2_ul.SetNpy(500)
    
    g2_xsecul.SetNpx(500)
    g2_xsecul.SetNpy(500)
    
    h2_ul_interp = g2_ul.GetHistogram().Clone()
    h2_ul_interp.SetName(f"{h2_ul.GetName()}_interp")
    h2_ul_interp.GetXaxis().SetTitle(h2_ul.GetXaxis().GetTitle())
    h2_ul_interp.GetYaxis().SetTitle(h2_ul.GetYaxis().GetTitle())
    h2_ul_interp.GetZaxis().SetTitle(h2_ul.GetZaxis().GetTitle())
    
    contours = get_contours_from_hist(
        hist = h2_ul_interp.Clone(),
        contour_value = 1.0,
        nsmooth = 0,
    )
    
    h2_xsecul_interp = g2_xsecul.GetHistogram().Clone()
    h2_xsecul_interp.SetName(f"{h2_xsecul.GetName()}_interp")
    h2_xsecul_interp.GetXaxis().SetTitle(h2_xsecul.GetXaxis().GetTitle())
    h2_xsecul_interp.GetYaxis().SetTitle(h2_xsecul.GetYaxis().GetTitle())
    h2_xsecul_interp.GetZaxis().SetTitle(h2_xsecul.GetZaxis().GetTitle())
    
    result = sortedcontainers.SortedDict()
    
    result["g2_ul"] = g2_ul
    result["h2_ul"] = h2_ul
    result["h2_ul_interp"] = h2_ul_interp
    
    result["g2_xsecul"] = g2_xsecul
    result["h2_xsecul"] = h2_xsecul
    result["h2_xsecul_interp"] = h2_xsecul_interp
    
    #result["d_xsecul_per_ctau"] = d_xsecul_per_ctau
    result["contour"] = None
    
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
        "--extratext",
        help = "Extra text to be printed on canvas",
        type = str,
        required = False,
    )
    
    parser.add_argument(
        "--xsecfile",
        help = "Cross section csv file",
        type = str,
        required = True,
    )
    
    parser.add_argument(
        "--cmsextratext",
        help = "Will set the CMS extra text; can be empty string \"\"",
        type = str,
        required = False,
        default = "Preliminary",
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
    
    
    arr_xsec_data = numpy.loadtxt(args.xsecfile, delimiter = ",")
    
    #nelements = arr_data.shape[0]
    #arr_mass = numpy.array(arr_data[:, 0], dtype = numpy.float64)
    #arr_xsec = numpy.array(arr_data[:, 1], dtype = numpy.float64)
    
    d_xsecs = {_row[0]: {
        "nom": _row[1],
        "unc_d": _row[1],
        "unc_u": _row[2]
    } for _row in arr_xsec_data}
    
    print(d_xsecs)
    
    d_limits = sortedcontainers.SortedDict()
    
    #for ipoint, injson in enumerate(args.jsons[: 20]) :
    for ipoint, injson in enumerate(args.jsons) :
        
        limits = cmut.load_config(injson)["120.0"]
        #print(limits)
        
        sig_info = cmut.parse_stau_samplestring(
            s = injson,
            regexp = "SMS-TStauStau_MStau-(?P<mstau>\w+)_ctau-(?P<ctau>\w+)_mLSP-(?P<mlsp>\d+)",
        )
        
        mstau_str = sig_info["mstau"]
        ctau_str = sig_info["ctau"]
        
        limits["mstau_str"] = mstau_str
        limits["ctau_str"] = ctau_str
        
        mstau = float(mstau_str)
        ctau = float(ctau_str[0: -2].replace("p", "."))
        
        d_limits[(mstau, ctau)] = limits
    
    
    outfile = ROOT.TFile.Open(args.output, "RECREATE")
    
    # Contours
    results = sortedcontainers.SortedDict()
    
    results["obs"] = get_exclusion(
        d_limits = d_limits,
        d_xsecs = d_xsecs,
        limit_key = "obs",
        name = "obs",
        axis_label = "Observed UL",
    )
    
    results["exp"] = get_exclusion(
        d_limits = d_limits,
        d_xsecs = d_xsecs,
        limit_key = "exp0",
        name = "exp",
        axis_label = "Expected UL at 95% CL",
    )
    
    results["exp_p1"] = get_exclusion(
        d_limits = d_limits,
        d_xsecs = d_xsecs,
        limit_key = "exp+1",
        name = "exp_p1",
        axis_label = "Expected UL (+1#sigma) at 95% CL",
    )
    
    results["exp_m1"] = get_exclusion(
        d_limits = d_limits,
        d_xsecs = d_xsecs,
        limit_key = "exp-1",
        name = "exp_m1",
        axis_label = "Expected UL (-1#sigma) at 95% CL",
    )
    
    #results["exp_p2"] = get_exclusion(
    #    d_limits = d_limits,
    #    d_xsecs = d_xsecs,
    #    limit_key = "exp+2",
    #    name = "exp_p2",
    #    axis_label = "Expected UL (+2#sigma) at 95% CL",
    #)
    #
    #results["exp_m2"] = get_exclusion(
    #    d_limits = d_limits,
    #    d_xsecs = d_xsecs,
    #    limit_key = "exp-2",
    #    name = "exp_m2",
    #    axis_label = "Expected UL (-2#sigma) at 95% CL",
    #)
    
    
    # 1D exclusions
    d_xsecul_per_ctau = sortedcontainers.SortedDict()
    
    for (mstau, ctau), limits in d_limits.items() :
        
        mstau_str = limits["mstau_str"]
        ctau_str = limits["ctau_str"]
        
        ctau_key = (ctau, ctau_str)
        print(mstau_str, ctau_str)
        
        l_limit_keys = ["obs", "exp0", "exp+1", "exp-1", "exp+2", "exp-2"]
        skip = False
        
        for limit_key in l_limit_keys :
            if limit_key not in limits :
                
                print(f"Warning: could not find limit {limit_key} for mstau={mstau}, ctau={ctau}. Skipping this signal point.")
                skip = True
                break
        
        if skip :
            continue
        
        if ctau_key not in d_xsecul_per_ctau :
            
            d_xsecul_per_ctau[ctau_key] = {
                "g1_xsec_theory": ROOT.TGraph(),
                "g1_xsecul_obs": ROOT.TGraph(),
                "g1_xsecul_exp": ROOT.TGraph(),
                "g1_xsecul_exp_pm1": ROOT.TGraphAsymmErrors(),
                "g1_xsecul_exp_pm2": ROOT.TGraphAsymmErrors(),
            }
            
            for key, val in d_xsecul_per_ctau[ctau_key].items() :
                
                val.SetName(f"{key}_{ctau_str}")
                val.SetTitle(f"{key}_{ctau_str}")
        
        xsec = d_xsecs[mstau]["nom"]
        ipoint = d_xsecul_per_ctau[ctau_key]["g1_xsecul_exp"].GetN()
        
        d_xsecul_per_ctau[ctau_key]["g1_xsec_theory"].SetPoint(ipoint, mstau, xsec)
        d_xsecul_per_ctau[ctau_key]["g1_xsecul_obs"].SetPoint(ipoint, mstau, limits["obs"]*xsec)
        d_xsecul_per_ctau[ctau_key]["g1_xsecul_exp"].SetPoint(ipoint, mstau, limits["exp0"]*xsec)
        
        d_xsecul_per_ctau[ctau_key]["g1_xsecul_exp_pm1"].SetPoint(ipoint, mstau, limits["exp0"]*xsec)
        d_xsecul_per_ctau[ctau_key]["g1_xsecul_exp_pm1"].SetPointError(
            ipoint,
            0.0, 0.0,
            abs(limits["exp-1"]-limits["exp0"])*xsec, abs(limits["exp+1"]-limits["exp0"])*xsec
        )
        
        d_xsecul_per_ctau[ctau_key]["g1_xsecul_exp_pm2"].SetPoint(ipoint, mstau, limits["exp0"]*xsec)
        d_xsecul_per_ctau[ctau_key]["g1_xsecul_exp_pm2"].SetPointError(
            ipoint,
            0.0, 0.0,
            abs(limits["exp-2"]-limits["exp0"])*xsec, abs(limits["exp+2"]-limits["exp0"])*xsec
        )
    
    
    print(f"Writing to: {args.output}")
    outfile.cd()
    
    for result_key, result in results.items() :
        
        outfile.cd()
        outfile.mkdir(result_key)
        outfile.cd(result_key)
        
        for key, val in result.items() :
            
            if (val is None) :
                
                continue
            
            print(f"Writing [{result_key}] [{key}]")
            val.Write()
    
    canvas_name = "limits-vs-ctau0-mstau"
    
    CMS.setCMSStyle()
    CMS.SetExtraText(args.cmsextratext)
    CMS.SetLumi("")
    CMS.SetEnergy("")
    CMS.SetCMSPalette()
    CMS.ResetAdditionalInfo()
    
    #NRGBs = 5
    #NCont = 255
    #l_stops = [0.00, 0.15, 0.27, 0.55, 1.00]
    #l_red   = [0.50, 0.50, 1.00, 1.00, 1.00]
    #l_green = [0.50, 1.00, 1.00, 0.60, 0.50]
    #l_blue  = [1.00, 1.00, 0.50, 0.40, 0.50]
    #
    #ROOT.TColor.CreateGradientColorTable(
    #    NRGBs,
    #    numpy.array(l_stops),
    #    numpy.array(l_red),
    #    numpy.array(l_green),
    #    numpy.array(l_blue),
    #    NCont
    #)
    #
    #ROOT.gStyle.SetNumberContours(NCont)
    
    #ROOT.gStyle.SetPadTickX(0)
    #ROOT.gStyle.SetPadTickY(0)
    
    ymin = 1
    ymax = max([
        max(numpy.array(results[_key]["contour"].GetY()))
        if results[_key]["contour"] else ymin
        for _key in ["exp", "exp_m1", "exp_p1"]
    ])
    #ymax = (round(ymax/100)+4)*100
    #ymax = (round(ymax/100)+4)*1000
    ymax = 10**(numpy.log10(ymax)+2)
    
    canvas = CMS.cmsCanvas(
        canvName = canvas_name,
        x_min = 90, x_max = 600,
        y_min = ymin, y_max = ymax,
        nameXaxis = "m_{#tilde{#tau}} [GeV]",
        nameYaxis = "c#tau_{0}(#tilde{#tau}) [mm]",
        #square = CMS.kSquare,
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
    
    h2_colz = results["exp"]["h2_xsecul_interp"]
    
    zmin = 10**-1
    #zmin = h2_colz.GetMinimum()
    #zmin = 10**(round(numpy.log10(zmin))-1)
    
    zmax = 10**3
    #zmax = h2_colz.GetMaximum()
    #zmax = 10**(round(numpy.log10(zmax))+1)
    
    h2_colz.GetZaxis().SetRangeUser(zmin, zmax)
    h2_colz.GetZaxis().SetTitle("95% CL upper limit on cross section [fb]")
    h2_colz.GetZaxis().SetTitleSize(0.05)
    h2_colz.GetZaxis().SetTitleOffset(1.2)
    h2_colz.GetZaxis().SetLabelSize(0.05)
    
    legend = CMS.cmsLeg(canvas.GetLeftMargin(), 0.65, 1-canvas.GetRightMargin(), 1-canvas.GetTopMargin(), textSize = 0.0425, columns = 2)
    #legend = ROOT.TLegend(canvas.GetLeftMargin(), 0.65, 1-canvas.GetRightMargin(), 1-canvas.GetTopMargin())
    legend.SetFillColor(0)
    legend.SetFillStyle(1001)
    legend.SetBorderSize(1)
    #legend.SetMargin(0.25)
    legend.SetTextAlign(22)
    
    CMS.cmsDraw(results["exp"]["h2_xsecul_interp"], "colz")
    CMS.cmsDraw(results["exp"]["contour"], "sameL", lstyle = ROOT.kSolid, lcolor = ROOT.kRed+1, lwidth = 3)
    legend.AddEntry(results["exp"]["contour"], "Expected", "L")
    
    if (results["exp_p1"]["contour"] and results["exp_m1"]["contour"]) :
        
        CMS.cmsDraw(results["exp_p1"]["contour"], "sameL", lstyle = ROOT.kDashed, lcolor = ROOT.kRed+1, lwidth = 3)
        CMS.cmsDraw(results["exp_m1"]["contour"], "sameL", lstyle = ROOT.kDashed, lcolor = ROOT.kRed+1, lwidth = 3)
        
        legend.AddEntry(results["exp_p1"]["contour"], "Expected #pm 1#sigma_{experiment}", "L")
    
    CMS.GetcmsCanvasHist(canvas).GetXaxis().SetNdivisions(6, 5, 0)
    
    canvas.SetLogy()
    
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
    
    #canvas.RedrawAxis()
    legend.Draw()
    canvas.Update()
    
    canvas.SetLogz()
    
    outfile.cd()
    canvas.Write()
    
    canvas_outfile = f"{outdir}/{canvas_name}"
    CMS.SaveCanvas(canvas, f"{canvas_outfile}.pdf")
    os.system(f"pdftoppm -cropbox -r 600 -png -singlefile {canvas_outfile}.pdf {canvas_outfile}")
    
    
    for (ctau, ctau_str), d_xsecul in d_xsecul_per_ctau.items() :
        
        outfile.cd()
        outfile.mkdir(ctau_str)
        outfile.cd(ctau_str)
        
        for key, val in d_xsecul.items() :
            
            print(f"Writing [{ctau_str}] [{key}]")
            #print(numpy.array(val.GetY(), dtype = numpy.float64))
            val.Write()
        
        g1_xsec_theory = d_xsecul["g1_xsec_theory"]
        g1_xsecul_obs = d_xsecul["g1_xsecul_obs"]
        g1_xsecul_exp = d_xsecul["g1_xsecul_exp"]
        g1_xsecul_exp_pm1 = d_xsecul["g1_xsecul_exp_pm1"]
        g1_xsecul_exp_pm2 = d_xsecul["g1_xsecul_exp_pm2"]
        
        xmin = 50
        xmax = 750
        
        #ymin = min([_ele.GetHistogram().GetMinimum() for _ele in d_xsecul.values()])
        ymin = g1_xsec_theory.GetHistogram().GetMinimum()
        ymin = 10**(round(numpy.log10(ymin))-1)
        #ymin = 10**-2
        
        #ymax = max([_ele.GetHistogram().GetMaximum() for _ele in d_xsecul.values()])
        ymax = g1_xsec_theory.GetHistogram().GetMaximum()
        ymax = 10**(round(numpy.log10(ymax))+2)
        
        canvas_name = f"limits-vs-mstau_{ctau_str}"
        
        CMS.setCMSStyle()
        CMS.SetExtraText(args.cmsextratext)
        CMS.SetLumi(137.62, round_lumi = True)
        CMS.SetEnergy("13")
        CMS.ResetAdditionalInfo()
        
        canvas = CMS.cmsCanvas(
            canvName = canvas_name,
            x_min = xmin, x_max = xmax,
            y_min = ymin, y_max = ymax,
            nameXaxis = "m_{#tilde{#tau}} [GeV]",
            nameYaxis = "Cross section [fb]",
            square = CMS.kSquare,
            iPos = 0,
            extraSpace = 0.03,
            with_z_axis = False,
            scaleLumi = None,
            yTitOffset = 1.25
        )
        
        CMS.cmsDraw(g1_xsecul_exp_pm2, "3L", fcolor = ROOT.TColor.GetColor("#85D1FBff"))
        CMS.cmsDraw(g1_xsecul_exp_pm1, "same3L", fcolor = ROOT.TColor.GetColor("#FFDF7Fff"))
        CMS.cmsDraw(g1_xsecul_exp, "sameL", lstyle = ROOT.kDashed, lcolor = ROOT.kBlack, lwidth = 3)
        CMS.cmsDraw(g1_xsec_theory, "sameL", lstyle = ROOT.kDotted, lcolor = ROOT.kRed, lwidth = 3)
        #CMS.cmsDraw(g1_xsecul_obs, "LP")
        
        g1_xsecul_exp_pm1.SetLineWidth(0)
        g1_xsecul_exp_pm2.SetLineWidth(0)
        
        legend = CMS.cmsLeg(0.5, 0.60, 0.95, 0.90, textSize = 0.04, columns = 1)
        #leg.AddEntry(g1_xsecul_obs, "Observed", "LP")
        legend.AddEntry(g1_xsec_theory, "Theory (NLO+NLL)", "L")
        legend.AddEntry(0, "#kern[-0.3]{#bf{95% CL upper limits}}", "")
        legend.AddEntry(g1_xsecul_exp, "Expected", "L")
        legend.AddEntry(g1_xsecul_exp_pm1, "Expected #pm 1#sigma_{experiment}", "F")
        legend.AddEntry(g1_xsecul_exp_pm2, "Expected #pm 2#sigma_{experiment}", "F")
        
        canvas.SetLogy()
        
        label = (
            "pp#rightarrow#tilde{#tau}#bar{#tilde{#tau}} , "
            "#tilde{#tau}#rightarrow#tau#tilde{G} , "
            "m_{#tilde{G}} = 1 GeV , "
            "c#tau_{0}(#tilde{#tau}) = " + ctau_str[0:-2].replace("p", ".") + " mm"
        )
        
        if args.extratext :
            
            label = f"#splitline{{{args.extratext}}}{{{label}}}"
        
        model_text = ROOT.TLatex()
        model_text.SetTextSize(0.04)
        model_text.DrawLatex(100, 5*ymin, label)
        
        canvas.Write()
        #canvas.SaveAs(f"{outdir}/{canvas_name}.pdf", close = False)
        #canvas.SaveAs(f"{outdir}/{canvas_name}.png")
        
        canvas_outfile = f"{outdir}/{canvas_name}"
        CMS.SaveCanvas(canvas, f"{canvas_outfile}.pdf")
        
        # ROOT fucks up when saving as png
        # Convert from the pdf instead
        os.system(f"pdftoppm -cropbox -r 600 -png -singlefile {canvas_outfile}.pdf {canvas_outfile}")
    
    return 0


if (__name__ == "__main__") :
    
    main()