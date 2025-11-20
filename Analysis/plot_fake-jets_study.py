#!/usr/bin/env python3

import os
import json
from hepunits import d
import particle

import ROOT
ROOT.gROOT.SetBatch(True)

from utils import cms_lumi as CMS_lumi
import utils.utils as utils
import utils.commonutils as cmut

def main() :
    
    #infilename = "output/fake-jets_study/fake-jets_study.root"
    #infilename = "output/fake-jets_study_w-charged-genparticles/fake-jets_study.root"
    #infilename = "output/fake-jets_study_w-charged-genparticles_w-bveto/fake-jets_study.root"
    infilename = "output/fake-jets_study_w-charged-genparticles_wo-genleptons_w-bveto/fake-jets_study.root"
    infile = ROOT.TFile.Open(infilename, "READ")
    
    outdir = os.path.dirname(infilename)
    os.system(f"mkdir -p {outdir}")
    
    d_processes = {}
    d_processes["QCD"] = [
        "QCD_Pt_30to50_TuneCP5_13TeV_pythia8",
        "QCD_Pt_50to80_TuneCP5_13TeV_pythia8",
        "QCD_Pt_80to120_TuneCP5_13TeV_pythia8",
        "QCD_Pt_120to170_TuneCP5_13TeV_pythia8",
        "QCD_Pt_170to300_TuneCP5_13TeV_pythia8",
        "QCD_Pt_300to470_TuneCP5_13TeV_pythia8",
        "QCD_Pt_470to600_TuneCP5_13TeV_pythia8",
        "QCD_Pt_600to800_TuneCP5_13TeV_pythia8",
        "QCD_Pt_800to1000_TuneCP5_13TeV_pythia8",
        "QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8",
        "QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8",
        "QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8",
        "QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8",
        "QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8",
    ]
    
    l_processes = list(d_processes.keys())
    
    d_proc_colors = {}
    d_proc_colors["QCD"] = ROOT.TColor.GetColor("#92dadd")
    
    d_proc_linestyles = {}
    d_proc_linestyles["QCD"] = 1
    
    d_proc_labels = {}
    d_proc_labels["QCD"] = "QCD"
    
    d_sample_xsecs = {}
    
    with open("configs/crosssections.json", "r") as fopen:
        d_sample_xsecs = json.load(fopen)
    
    histname_nevents = "nEvents"
    
    l_histname = [
        "Jet_flavor",
        "Jet_dxy_vs_flavor",
    ]
    
    d_xrange = {}
    d_xrange["Jet_flavor"] = (1, 600)
    d_xrange["Jet_dxy_vs_flavor"] = (1, 600)
    
    d_yrange = {}
    d_yrange["Jet_flavor"] = (1e-6, 10)
    d_yrange["Jet_dxy_vs_flavor"] = (0.1, 100)
    
    d_xtitle = {}
    d_xtitle["Jet_flavor"] = "PDG ID"
    d_xtitle["Jet_dxy_vs_flavor"] = "PDG ID"
    
    d_ytitle = {}
    d_ytitle["Jet_flavor"] = "Fraction of jets"
    d_ytitle["Jet_dxy_vs_flavor"] = "Jet d_{xy} [cm]"
    
    d_hist = {}
    
    rebin = None
    #rebin = cmut.get_binning(0, 1, 0.04)
    
    latex = ROOT.TLatex()
    l_texts = []
    
    for proc, l_samples in d_processes.items() :
        
        d_hist[proc] = {}
        
        for histname in l_histname :
            
            d_hist[proc][histname] = None
            
            for sample in l_samples :
                print(f"[{sample}] [{histname}]")
                
                hist = cmut.get_hist(
                    histfile = infile,
                    histname = histname,
                    samples = [sample],
                    scales = [],
                    rebin = rebin,
                    underflow = True,
                    overflow = True,
                    set_min = None,
                    set_max = None,
                )
                
                hist_nevents = cmut.get_hist(
                    histfile = infile,
                    histname = histname_nevents,
                    samples = [sample],
                    scales = [],
                    rebin = rebin,
                    underflow = True,
                    overflow = True,
                    set_min = None,
                    set_max = None,
                )
                
                hist.Scale(d_sample_xsecs[sample] / hist_nevents.GetEntries())
                
                if d_hist[proc][histname] is None :
                    d_hist[proc][histname] = hist
                else :
                    d_hist[proc][histname].Add(hist)
                
                bin_start = d_hist[proc][histname].GetXaxis().FindBin(d_xrange[histname][0])
                bin_end = d_hist[proc][histname].GetXaxis().FindBin(d_xrange[histname][1])
            
            if isinstance(d_hist[proc][histname], ROOT.TH1D) :
                
                d_hist[proc][histname].Scale(1.0 / d_hist[proc][histname].Integral(bin_start, bin_end))
                
                for ibin in range(1, d_hist[proc][histname].GetNbinsX()+1) :
                    
                    bin_val = d_hist[proc][histname].GetBinContent(ibin)
                    pdgid = d_hist[proc][histname].GetBinLowEdge(ibin)
                    
                    latex_symbol = ""
                    
                    try :
                        latex_symbol = particle.Particle.from_pdgid(pdgid).latex_name.replace("\\", "#")
                    except Exception as excp:
                        pass
                    
                    if bin_val > d_yrange[histname][0] and pdgid :
                        
                        #text = ROOT.TText()
                        text = ROOT.TLatex()
                        #text.SetNDC(True)
                        #text.SetText(pdgid, 1.5*bin_val, f"{int(pdgid)}")
                        text.SetText(pdgid, 1.5*bin_val, latex_symbol)
                        text.SetTextAngle(45)
                        text.SetTextSize(0.025)
                        #text.SetTextFont(42)
                        text.SetTextColor(1)
                        l_texts.append(text)
            
            elif isinstance(d_hist[proc][histname], ROOT.TH2D) :
                
                bin_start_y = d_hist[proc][histname].GetYaxis().FindBin(d_yrange[histname][0])
                bin_end_y = d_hist[proc][histname].GetYaxis().FindBin(d_yrange[histname][1])
                
                d_hist[proc][histname].Scale(1.0 / d_hist[proc][histname].Integral(bin_start, bin_end, bin_start_y, bin_end_y))
    
    print(l_texts)
    
    for histname in l_histname :
        
        xrange = d_xrange[histname]
        yrange = d_yrange[histname]
        
        if isinstance(d_hist[l_processes[0]][histname], ROOT.TH1D) :
            
            l_hist = []
            outfilename = f"{outdir}/{histname}.pdf"
            
            for iproc, proc in enumerate(d_hist.keys()) :
                
                #color = utils.get_cms_colors(iproc)
                color = d_proc_colors[proc]
                linestyle = d_proc_linestyles[proc]
                
                print(proc, histname)
                hist = d_hist[proc][histname].Clone()
                hist.SetTitle(d_proc_labels[proc])
                hist.SetMarkerColor(color)
                hist.SetLineColor(color)
                hist.SetLineStyle(linestyle)
                hist.SetLineWidth(3)
                #hist.SetMarkerSize(0)
                hist.SetFillStyle(0)
                hist.SetDrawOption("hist")
                hist.SetOption("hist")
                
                l_hist.append(hist)
            
            utils.root_plot1D_legacy(
                l_hist = l_hist,
                l_text_overlay = l_texts,
                #ratio_num_den_pairs = [(_hist, l_hist[-1]) for _hist in l_hist[:-1]],
                outfile = outfilename,
                xrange = xrange,
                yrange = yrange,
                logx = False,
                logy = True,
                xtitle = d_xtitle[histname],
                ytitle = d_ytitle[histname],
                xtitle_ratio = d_xtitle[histname],
                ytitle_ratio = f"Ratio wrt {l_hist[-1].GetTitle()}",
                centertitlex = True, centertitley = True,
                centerlabelx = False, centerlabely = False,
                #gridx = True, gridy = True,
                ndivisionsx = None,
                stackdrawopt = "nostack",
                ratiodrawopt = "E1",
                ndivisionsy_ratio = (2, 5, 0),
                legendpos = "UL",
                #legendtitle = "[m_{#tilde{#tau}}(250), m_{LSP}(1)]",
                legendncol = 1,
                legendtextsize = 0.05,
                legendwidthscale = 1.5,
                #legendheightscale = 0.9,
                legendpadleft_extra = 0.05,
                CMSextraText = "Simulation Preliminary",
                #lumiText = "2018 (13 TeV)"
                lumiText = "13 TeV"
            )
        
        elif isinstance(d_hist[l_processes[0]][histname], ROOT.TH2D) :
            
            for iproc, proc in enumerate(d_hist.keys()) :
                
                outfilename = f"{outdir}/{histname}_{proc}.pdf"
                
                hist = d_hist[proc][histname].Clone()
                
                #ROOT.gROOT.LoadMacro("utils/tdrstyle.C")
                #ROOT.gROOT.ProcessLine("setTDRStyle()")
                #
                #ROOT.gROOT.SetStyle("tdrStyle")
                #ROOT.gROOT.ForceStyle(True)
                
                #canvas = ROOT.TCanvas("canvas", "canvas", 1000, 875)
                #canvas.UseCurrentStyle()
                
                canvas = utils.get_canvas()
                
                canvas.SetLeftMargin(0.13)
                canvas.SetRightMargin(0.225)
                canvas.SetTopMargin(0.08)
                canvas.SetBottomMargin(0.13)
                
                hist.GetXaxis().SetRangeUser(*xrange)
                hist.GetXaxis().SetTitle(d_xtitle[histname])
                hist.GetXaxis().SetTitleSize(0.05)
                hist.GetXaxis().SetTitleOffset(1.1)
                hist.GetXaxis().CenterTitle(True)
                hist.GetXaxis().SetLabelSize(0.045)
                #if (d_config["nDivX"] is not None) : hist.GetXaxis().SetNdivisions(*d_config["nDivX"], True)
                
                hist.GetYaxis().SetRangeUser(*yrange)
                hist.GetYaxis().SetTitle(d_ytitle[histname])
                hist.GetYaxis().SetTitleSize(0.05)
                hist.GetYaxis().SetTitleOffset(1.1)
                hist.GetYaxis().CenterTitle(True)
                hist.GetYaxis().SetLabelSize(0.045)
                #if (d_config["nDivY"] is not None) : hist.GetYaxis().SetNdivisions(*d_config["nDivY"], True)
                
                hist.GetZaxis().SetTitle("a.u.")
                hist.GetZaxis().SetTitleSize(0.05)
                hist.GetZaxis().SetTitleOffset(1.55)
                hist.GetZaxis().CenterTitle(True)
                hist.GetZaxis().SetLabelSize(0.045)
                
                hist.SetMinimum(1e-6)
                hist.SetMaximum(1)
                
                #utils.cpalette_nipy_spectral.set()
                
                hist.Draw("colz")
                
                canvas.SetLogy(True)
                canvas.SetLogz(True)
                
                CMS_lumi.CMS_lumi(pad = canvas, iPeriod = 0, iPosX = 0, CMSextraText = "Simulation Preliminary", lumiText = "13 TeV")
                
                outfile_noext, _ = os.path.splitext(outfilename)
                
                canvas.SaveAs(outfilename)
                canvas.Close()
                
                utils.logger.info(f"Converting: {[outfilename]} --> [{outfile_noext}.png]")
                os.system(f"pdftoppm -png -r 600 -cropbox -singlefile {outfilename} {outfile_noext}")


if __name__ == "__main__" :
    
    main()