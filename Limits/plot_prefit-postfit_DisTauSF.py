#!/usr/bin/env -S python3 -u

import argparse
import numpy
import os
import ROOT

import utils.commonutils as cmut
import utils.utils as utils

#>> Doing postfit: bin1_mumu_2018,Other
#>> Doing postfit: bin1_mumu_2018,ST
#>> Doing postfit: bin1_mumu_2018,TT
#>> Doing postfit: bin1_mumu_2018,TTX
#>> Doing postfit: bin1_mumu_2018,WLNu
#>> Doing postfit: bin1_mumu_2018,ZLL
#>> Doing postfit: bin1_mumu_2018,TotalBkg
#>> Doing postfit: bin1_mumu_2018,TotalSig
#>> Doing postfit: bin1_mumu_2018,TotalProcs

#>> Doing postfit: bin1_mutau_fail_2018,Other
#>> Doing postfit: bin1_mutau_fail_2018,QCD
#>> Doing postfit: bin1_mutau_fail_2018,ST
#>> Doing postfit: bin1_mutau_fail_2018,TT
#>> Doing postfit: bin1_mutau_fail_2018,TTX
#>> Doing postfit: bin1_mutau_fail_2018,WLNu
#>> Doing postfit: bin1_mutau_fail_2018,ZLL
#>> Doing postfit: bin1_mutau_fail_2018,sig
#>> Doing postfit: bin1_mutau_fail_2018,TotalBkg
#>> Doing postfit: bin1_mutau_fail_2018,TotalSig
#>> Doing postfit: bin1_mutau_fail_2018,TotalProcs

#bin1_mutau_fail_2018_prefit
#bin1_mutau_pass_2018_prefit
#bin2_mutau_fail_2018_prefit
#bin2_mutau_pass_2018_prefit
#bin1_mumu_2018_prefit

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
        "--era",
        help = "Era",
        type = str,
        required = True,
        choices = ["2016_preVFP", "2016_postVFP", "2017", "2018"],
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
        "--channels",
        help = "Channels",
        type = str,
        required = True,
        nargs = "+",
        choices = ["mumu", "mutau"],
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    #infile_name = "/home/soham/nfs_dust/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/tmp/test_DisTauSF_mass/DisTauSF/channels_all/eras_all/ZMT_wp-p80_dxy-gt-0p00/postfit_s.root"
    #infile_name = "/home/soham/nfs_dust/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/tmp/test_DisTauSF_mass/DisTauSF/channels_all/eras_all/ZMT_wp-p99_dxy-gt-0p00/postfit_s.root"
    #infile_name = f"/home/soham/nfs_dust/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/tmp/test_DisTauSF_mass/DisTauSF/channels_all/eras_all/{dir_name}/postfit_s.root"
    
    d_channel_info = {}
    
    d_channel_info["mumu"] = {
        "label": "#mu#mu",
        "nbins": 1,
        #"eras": ["2018"],
        "eras": [args.era],
        "categories": [""],
        "data": "data_obs",
        "processes_bkg": ["ZLL", "WLNu", "TTX", "TT", "ST", "Other"],
        "processes_sig": [],
        "bkg_total": "TotalBkg",
        "sig_total": "TotalSig",
        "procs_total": "TotalProcs",
        "fits": ["prefit", "postfit"],
    }
    
    d_channel_info["mutau"] = {
        "label": "#mu#tau_{h}",
        "nbins": 1,
        #"eras": ["2018"],
        "eras": [args.era],
        "categories": ["pass", "fail"],
        "data": "data_obs",
        #"processes_bkg": ["ZLL", "WLNu", "TTX", "TT", "ST", "QCD", "Other"],
        "processes_bkg": ["Other", "QCD", "WLNu", "Top", "ZLL"],
        "processes_sig": ["sig"],
        "bkg_total": "TotalBkg",
        "sig_total": "TotalSig",
        "procs_total": "TotalProcs",
        "fits": ["prefit", "postfit"],
    }
    
    #["#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6", "#a96b59", "#e76300", "#b9ac70", "#717581", "#92dadd"]
    d_proc_info = {
        "sig": {
            "label": "DY (#mu#tau_{h})",
            "color": "#3f90da",
        },
        "ZLL": {
            "label": "DY (other)",
            "color": "#ffa90e",
        },
        "WLNu": {
            "label": "W+jets",
            #"color": "#bd1f01",
            "color": "#b9ac70",
        },
        #"TTX": {
        #    "label": "t#bar{t}X",
        #    "color": "#94a4a2",
        #},
        #"TT": {
        #    "label": "t#bar{t}",
        #    "color": "#832db6",
        #},
        #"ST": {
        #    "label": "single t",
        #    "color": "#a96b59",
        #},
        "Top": {
            "label": "Top",
            "color": "#a96b59",
        },
        "QCD": {
            "label": "Multijet",
            "color": "#e76300",
        },
        "Other": {
            "label": "Other SM",
            #"color": "#b9ac70",
            "color": "#94a4a2",
        },
    }
    
    for infile_name in args.input :
        
        outdir = args.outdir if (args.outdir) else os.path.dirname(infile_name)
        os.system(f"mkdir -p {outdir}")
        
        parsed_result = cmut.parse_string_regex(
            s = infile_name,
            regexp = "ZMT_wp-(?P<wp>\w+)_dxy-gt-(?P<dxy>\w+)",
        )
        
        wp = parsed_result["wp"]
        dxy = parsed_result["dxy"]
        
        wp_key = f"wp-{wp}"
        dxy_key = f"dxy-gt-{dxy}"
        
        wp = float(wp.replace("p", "."))
        dxy = float(dxy.replace("p", "."))
        
        sig_with_bkg = True
        
        print(f"Opening file: {infile_name}")
        infile = ROOT.TFile.Open(infile_name)
        
        #for ch_name, ch_info in d_channel_info.items() :
        for ch_name in args.channels :
            
            ch_info = d_channel_info[ch_name]
            
            for fit in ch_info["fits"] :
                
                for era in ch_info["eras"] :
                    
                    for cat in ch_info["categories"] :
                        
                        d_hist = {}
                        
                        procs = ch_info["processes_bkg"] + ch_info["processes_sig"] + [
                            ch_info["data"],
                            ch_info["bkg_total"],
                            ch_info["sig_total"],
                            ch_info["procs_total"],
                        ]
                        
                        for iproc, proc in enumerate(procs) :
                            
                            nbins = ch_info["nbins"]
                            hist_proc_name = f"{ch_name}_{era}_{cat}_{proc}_{fit}"
                            hist_proc_name = hist_proc_name.replace("__", "_")
                            print(f"Getting histogram: {hist_proc_name}")
                            hist_proc = ROOT.TH1F(hist_proc_name, proc, nbins, 0, nbins)
                            hist_proc.SetDirectory(0)
                            
                            #color = utils.ColorIterator(iproc, 0)
                            color = 1
                            label = hist_proc.GetTitle()
                            
                            if (proc in d_proc_info) :
                                
                                color = ROOT.TColor.GetColor(d_proc_info[proc]["color"])
                                label = d_proc_info[proc]["label"]
                            
                            for ibin in range(nbins) :
                                
                                inhist_name = f"bin{ibin+1}_{ch_name}_{cat}_{era}_{fit}/{proc}"
                                inhist_name = inhist_name.replace("__", "_")
                                print(inhist_name)
                                inhist = infile.Get(inhist_name).Clone()
                                inhist.SetDirectory(0)
                                
                                bin_val = inhist.GetBinContent(1)
                                bin_err = inhist.GetBinError(1)
                                
                                hist_proc.SetBinContent(ibin+1, bin_val)
                                hist_proc.SetBinError(ibin+1, bin_err)
                                
                            
                            hist_proc.SetTitle(label)
                            hist_proc.SetOption("hist")
                            hist_proc.SetFillColor(color)
                            hist_proc.SetLineColor(color)
                            hist_proc.SetMarkerColor(color)
                            hist_proc.SetLineWidth(0)
                            #hist_proc.SetMarkerSize(0)
                            
                            print(hist_proc_name)
                            hist_proc.Print()
                            
                            d_hist[proc] = hist_proc#.Clone()
                        
                        hist_data = d_hist[ch_info["data"]]
                        hist_data.SetTitle("Data")
                        hist_data.SetLineColor(1)
                        hist_data.SetMarkerColor(1)
                        hist_data.SetMarkerSize(2)
                        hist_data.SetMarkerStyle(20)
                        hist_data.SetLineWidth(2)
                        hist_data.SetOption("PE1")
                        hist_data.SetFillStyle(0)
                        hist_data.SetFillColor(0)
                        
                        l_hist_bkg = [d_hist[_proc] for _proc in ch_info["processes_bkg"]]
                        l_hist_sig = [d_hist[_proc] for _proc in ch_info["processes_sig"]]
                        
                        l_hist = l_hist_bkg
                        hist_procs_total = d_hist[ch_info["bkg_total"]]
                        
                        if (sig_with_bkg) :
                            
                            l_hist = l_hist_bkg + l_hist_sig
                            
                            hist_procs_total = d_hist[ch_info["procs_total"]]
                        
                        hist_procs_total.SetFillColorAlpha(1, 1.0)
                        hist_procs_total.SetFillStyle(3354)
                        hist_procs_total.SetMarkerSize(0)
                        hist_procs_total.SetMarkerStyle(0)
                        hist_procs_total.SetLineWidth(0)
                        #hist_procs_total.SetLineColor(2)
                        hist_procs_total.SetOption("E2")
                        
                        hist_procs_total.SetTitle("Pred. uncertainty")
                        
                        bin_max = max(hist_data.GetBinContent(hist_data.GetMaximumBin()), hist_procs_total.GetBinContent(hist_procs_total.GetMaximumBin()))
                        
                        yrange_power = int(numpy.log10(bin_max))+1
                        
                        title = f"WP={wp}, d_{{xy}}>{dxy} cm, {ch_info['label']}, {cat}, {fit}"
                        
                        plotfile = f"{outdir}/{wp_key}_{dxy_key}_{ch_name}_{era}_{cat}_{fit}.pdf".replace("__", "_")
                        
                        utils.root_plot1D_legacy(
                            l_hist = l_hist,
                            l_hist_overlay = [hist_procs_total, hist_data],
                            #l_hist_overlay = [hist_data],
                            #outfile = f"{indir}/{ch_name}_{era}_{cat}_{fit}.pdf".replace("__", "_"),
                            outfile = plotfile,
                            xrange = (hist_data.GetXaxis().GetXmin(), hist_data.GetXaxis().GetXmax()),
                            #yrange = (10**(yrange_power-3), 10**(yrange_power+5)),
                            yrange = (1e-1, 10**(yrange_power+5)),
                            #yrange = (1e-3, 1e10),
                            ratio_num_den_pairs = [(hist_data, hist_procs_total)],
                            ratio_mode = "data",
                            #no_xerror = True,
                            logx = False, logy = True,
                            title = "",
                            xtitle = "", ytitle = "Events",
                            xtitle_ratio = "", ytitle_ratio = "Data/Total pred.",
                            yrange_ratio = (0, 2),
                            centertitlex = True, centertitley = True,
                            centerlabelx = False, centerlabely = False,
                            gridx = False, gridy = False,
                            ndivisionsx = (ch_info["nbins"], 1, 0),
                            ndivisionsy = None,
                            ndivisionsy_ratio = (4, 5, 0),
                            stackdrawopt = "",
                            ratiodrawopt = "E1",
                            legendpos = "UR",
                            legendncol = 2,
                            legendtextsize = 0.045,
                            legendtitle = title,
                            legendheightscale = 0.65, legendwidthscale = 1.5,
                            #CMSextraText = "Preliminary",
                            CMSextraText = "Private Work",
                            lumiText = f"{args.era} (13 TeV)"
                        )
                        
                        #os.system(f"pdftoppm -png -r 600 -cropbox -singlefile {plotfile} {plotfile}")
        
        infile.Close()
    
    return 0


if __name__ == "__main__" :
    
    main()