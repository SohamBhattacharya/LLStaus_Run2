#!/usr/bin/env python3

import argparse
import functools
import numpy
import operator

import ROOT
import CombineHarvester.CombineTools.ch as ch

import utils.commonutils as cmut


def set_stat_errs(
    ch_obj,
    target,
    hist,
    bin_ids,
    bin_names,
) :
    
    for bin_id, bin_name in zip(bin_ids, bin_names) :
        
        #print(bin_id, bin_name)
        val = hist.GetBinContent(bin_id)
        val = 0 if (val < 0) else val
        err = hist.GetBinError(bin_id)
        err_rel = numpy.clip(1.0+(err/val), 0.001, 2.0) if val else 2.0
        
        ch_obj.cp().bin([bin_name]).AddSyst(
            target = target,
            #name = "stat_$PROCESS_$BIN_$CHANNEL_$ERA",
            name = "stat_$PROCESS_$BIN",
            type = "lnN",
            valmap = ch.SystMap()(
                (err_rel)
            )
        )


def set_syst_errs(
    ch_obj,
    target,
    hist_nom,
    hist_u,
    hist_d,
    bin_ids,
    bin_names,
    systname,
    iseracorr,
) :
    
    for bin_id, bin_name in zip(bin_ids, bin_names) :
        
        val = hist_nom.GetBinContent(bin_id)
        val = 0 if (val < 0) else val
        err_u = hist_u.GetBinContent(bin_id) - val
        err_d = hist_d.GetBinContent(bin_id) - val
        
        #err_u_rel = numpy.clip(1.0+(err_u/val), 0.001, 2.0) if val else 2.0
        #err_d_rel = numpy.clip(1.0+(err_d/val), 0.001, 2.0) if val else 0.001
        
        err_u_rel = numpy.clip(1.0+(err_u/val), 0.001, 2.0) if val else 1.01
        err_d_rel = numpy.clip(1.0+(err_d/val), 0.001, 2.0) if val else 0.99
        
        name_tmp = f"syst_{systname}"
        
        if (not iseracorr) :
            name_tmp = f"{name_tmp}_$ERA"
        
        ch_obj.cp().bin([bin_name]).AddSyst(
            target = target,
            name = name_tmp,
            type = "lnN",
            valmap = ch.SystMap()(
                (err_d_rel, err_u_rel)
            )
        )


def get_and_set_syst_errs(
    procname,
    inhistfile,
    systs,
    samples,
    scales,
    rebin,
    bin_ids,
    bin_names,
    hist_nom,
    ch_obj,
    target,
) :
    
    d_hist_proc_syst = {}
    
    for systname, systinfo in systs.items() :
        
        d_hist_proc_syst[systname] = {}
        
        for systvar in ["u", "d"] :
            
            if (isinstance(systinfo[systvar], str)) :
                
                hist_proc_syst = cmut.get_hist(
                    histfile = inhistfile,
                    histname = systinfo[systvar],
                    samples = samples,
                    scales = scales,
                    rebin = rebin,
                    set_min = 0,
                )
            
            elif (isinstance(systinfo[systvar], dict)) :
                print(systname, systinfo[systvar], procname)
                
                histname = None
                
                if (procname in systinfo[systvar]) :
                    
                    histname = systinfo[systvar][procname]
                
                elif ("*" in systinfo[systvar]) :
                    
                    histname = systinfo[systvar]["*"]
                
                if not histname :
                    
                    cmut.logger.error(
                        f"Could not find process {procname} in systematics {systname}: {systvar}: {systinfo[systvar]}. "
                    )
                
                hist_proc_syst = cmut.get_hist(
                    histfile = inhistfile,
                    histname = histname,
                    samples = samples,
                    scales = scales,
                    rebin = rebin,
                    #set_min = 0,
                    set_min = 0.001,
                )
                
            elif (isinstance(systinfo[systvar], float)) :
                
                hist_proc_syst = hist_nom.Clone()
                hist_proc_syst.Scale(systinfo[systvar])
            
            else :
                
                cmut.logger.error(
                    f"Could not decipher systematics {systname}: {systvar}: {systinfo[systvar]}. "
                    "Must be one of [str(histname), dict of {procname: str(histname)}, float]."
                )
            
            #hist_proc_syst.Print("range")
            d_hist_proc_syst[systname][systvar] = hist_proc_syst
        
        set_syst_errs(
            ch_obj = ch_obj,
            target = target,
            hist_nom = hist_nom,
            hist_u = d_hist_proc_syst[systname]["u"],
            hist_d = d_hist_proc_syst[systname]["d"],
            bin_ids = bin_ids,
            bin_names = bin_names,
            systname = systname,
            iseracorr = systinfo["iseracorr"],
        )


def set_rateparams(
    rateparams,
    #channel,
    #procname,
    ch_obj,
    target,
) :
    
    for paramname, paraminfo in rateparams.items() :
        
        name = f"rp_{paramname}"
        
        if (not paraminfo["iseracorr"]) :
            
            name = f"{name}_$ERA"
        
        ch_obj.cp().AddSyst(
            target = target,
            name = name,
            type = "rateParam",
            valmap = ch.SystMap()(
                (paraminfo["init"])
                #(paraminfo["init"], paraminfo["min"], paraminfo["max"])
            )
        )


def main() :
    
    # Argument parser
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        "--configs",
        help = "Configuration yaml files (per channel, per era). A (channel, era) combination must not be repeated.",
        type = str,
        nargs = "+",
        required = True,
    )
    
    parser.add_argument(
        "--outdir",
        help = "Output directory",
        type = str,
        required = False,
        default = "tmp/test_cards/combine"
    )
    
    parser.add_argument(
        "--chcombos",
        help = "Custom channel combinations to create cards for; comma separated lists: ch1,ch2 ch3,ch5",
        type = str,
        nargs = "*",
        required = False,
        default = [],
    )
    
    # Parse arguments
    args = parser.parse_args()
    #d_args = vars(args)
    
    
    d_config_all = {}
    cb = ch.CombineHarvester()
    
    for cfg in args.configs :
        
        d_config = cmut.load_config(cfg)
        #print(d_config)
        
        era = d_config["era"]
        channel = d_config["channel"]
        cfg_key = (channel, era)
        
        assert(cfg_key not in d_config_all)
        
        d_config_all[cfg_key] = d_config
    
    hist_template = None
    
    for cfg_key, d_config in d_config_all.items() :
        
        print("*"*50)
        print(cfg_key)
        print("*"*50)
        
        rebin = None
        nbins = d_config["nbins"]
        
        if (d_config["rebin"] is not None) :
            
            rebin = numpy.array(d_config["rebin"], dtype = numpy.float64)
            
            if (nbins <= 0) :
                
                nbins = len(rebin)-1
        
        assert(nbins > 0)
        
        analysis = d_config["analysis"]
        era = d_config["era"]
        lumi = eval(d_config["lumi"]) if isinstance(d_config["lumi"], str) else d_config["lumi"]
        channel = d_config["channel"]
        binstart = d_config["binstart"]
        categories = [(_binnum, f"bin{_binnum}_{channel}_{era}") for _binnum in range(binstart, binstart+nbins)]
        bin_ids = [_cat[0] for _cat in categories]
        bin_names = [_cat[1] for _cat in categories]
        print(bin_ids, bin_names)
        
        if ("sig" in d_config) :
            
            if ("loadconfig" in d_config["sig"]) :
                
                d_config["sig"].update(cmut.load_config(d_config["sig"]["loadconfig"]))
            
            cb.AddProcesses(
                mass = list(d_config["sig"]["procs"].keys()),
                analysis = [analysis],
                era = [era],
                channel = [channel],
                procs = ["sig"],
                bin = categories,
                signal = True,
            )
        
        cb.AddProcesses(
            mass = ["*"],
            analysis = [analysis],
            era = [era],
            channel = [channel],
            procs = list(d_config["bkg"]["procs"].keys()),
            bin = categories,
            signal = False,
        )
        
        if ("obs" in d_config) :
            
            cb.AddObservations(
                mass = ["*"],
                analysis = [analysis],
                era = [era],
                channel = [channel],
                bin = categories,
            )
            
            samples = d_config["obs"]["samples"]
            scale = eval(d_config["obs"]["scaleby"])
            inhistfile_obs = ROOT.TFile.Open(d_config["obs"]["histfile"])
            
            hist_obs = cmut.get_hist(
                histfile = inhistfile_obs,
                histname = d_config["obs"]["nominal"],
                samples = samples,
                rebin = rebin,
                set_min = 0,
            )
            
            cb.cp().channel([channel]).era([era]).ForEachObs(
                lambda x : x.set_rate(
                    scale * hist_obs.GetBinContent(x.bin_id())
            ))
        
        d_xsec = cmut.load_config(d_config["bkg"]["xsecfile"])
        #print(d_xsec)
        
        if ("sig" in d_config) :
            
            inhistfile_sig = ROOT.TFile.Open(d_config["sig"]["histfile"])
            cutflows_sig = cmut.load_config(d_config["sig"]["cutflowsfile"])
            
            for procname, procinfo in d_config["sig"]["procs"].items() :
                
                samples = procinfo["samples"]
                neventkey = d_config["sig"]["neventkey"]
                issusy = procinfo["issusy"]
                ismc = True if issusy else procinfo["ismc"]
                
                # Parse the signal sample name to the the stau mass
                # Get the corresponding xsec
                # Update the xsec dictionary
                if (issusy) :
                    
                    d_xsec.update(cmut.get_stau_xsec_dict(
                        l_samplestr = samples,
                        xsecfile = d_config["sig"]["xsecfile"],
                        regexp = "SMS-TStauStau_MStau-(?P<mstau>\d+)_ctau-(?P<ctau>\w+)_mLSP-(?P<mlsp>\d+)",
                    ))
                
                scales = []
                
                for sample in samples :
                    
                    scale = eval(procinfo["scaleby"])
                    
                    if (ismc) :
                        
                        scale *= lumi
                        
                        if (procinfo["xsnorm"]) :
                            
                            neventtot = functools.reduce(operator.getitem, [sample]+neventkey.split("."), cutflows_sig)
                            scale *= (d_xsec[sample] / neventtot)
                    
                    scales.append(scale)
                
                #print(scales)
                
                hist_proc_nom = cmut.get_hist(
                    histfile = inhistfile_sig,
                    histname = d_config["sig"]["nominal"],
                    samples = samples,
                    scales = scales,
                    rebin = rebin,
                    set_min = 0,
                )
                
                #hist_proc_nom.Print("range")
                
                cb.cp().mass([procname]).channel([channel]).process(["sig"]).era([era]).ForEachProc(
                    lambda x : x.set_rate(
                        hist_proc_nom.GetBinContent(x.bin_id()) if hist_proc_nom.GetBinContent(x.bin_id()) else 0.001
                ))
                
                set_stat_errs(
                    ch_obj = cb.cp().mass([procname]).process(["sig"]).era([era]),
                    target = cb,
                    hist = hist_proc_nom,
                    bin_ids = bin_ids,
                    bin_names = bin_names,
                )
                
                procinfo["systematics"] = {}
                
                for systkey, systinfo in d_config["systematics"].items() :
                    
                    if (procname in systinfo["procs"]["sig"] or "sig" in systinfo["procs"]["sig"]) :
                        
                        procinfo["systematics"][systkey] = systinfo
                
                if (len(procinfo["systematics"])) :
                    
                    get_and_set_syst_errs(
                        procname = procname,
                        inhistfile = inhistfile_sig,
                        systs = procinfo["systematics"],
                        samples = samples,
                        scales = scales,
                        rebin = rebin,
                        bin_ids = bin_ids,
                        bin_names = bin_names,
                        hist_nom = hist_proc_nom,
                        ch_obj = cb.cp().mass([procname]).process(["sig"]).era([era]),
                        target = cb,
                    )
                
                procinfo["rateparams"] = {}
                
                for rpkey, rpinfo in d_config.get("rateparams", {}).items() :
                    
                    if (procname in rpinfo["procs"]["sig"]) :
                        
                        procinfo["rateparams"][rpkey] = rpinfo
                
                if (len(procinfo["rateparams"])) :
                    
                    set_rateparams(
                        #channel = channel,
                        #procname = "sig",
                        rateparams = procinfo["rateparams"],
                        #ch_obj = cb.cp().mass([procname]).process(["sig"]).era([era]),
                        ch_obj = cb.cp().mass([procname]).era([era]),
                        target = cb,
                    )
            
            inhistfile_sig.Close()
        
        
        for procname, procinfo in d_config["bkg"]["procs"].items() :
            
            ismc = procinfo["ismc"]
            #scale = eval(procinfo["scaleby"])
            inhistfile = ROOT.TFile.Open(procinfo["histfile"])
            cutflows_bkg = cmut.load_config(d_config["bkg"]["cutflowsfile"])
            samples = procinfo["samples"]
            #scales = [scale]*len(samples)
            
            scales = []
            
            for sample in samples :
                
                scale = eval(procinfo["scaleby"])
                
                if (ismc) :
                    
                    scale *= lumi
                    
                    if (procinfo["xsnorm"]) :
                        
                        neventtot = functools.reduce(operator.getitem, [sample]+neventkey.split("."), cutflows_bkg)
                        scale *= (d_xsec[sample] / neventtot)
                
                scales.append(scale)
            
            hist_proc_nom = cmut.get_hist(
                histfile = inhistfile,
                histname = procinfo["nominal"],
                samples = samples,
                scales = scales,
                rebin = rebin,
                set_min = 0,
            )
            
            #hist_proc_nom.Print("range")
            
            cb.cp().channel([channel]).process([procname]).era([era]).ForEachProc(
                lambda x : x.set_rate(
                    hist_proc_nom.GetBinContent(x.bin_id()) if hist_proc_nom.GetBinContent(x.bin_id()) else 0.001
            ))
            
            set_stat_errs(
                ch_obj = cb.cp().process([procname]).era([era]),
                target = cb,
                hist = hist_proc_nom,
                bin_ids = bin_ids,
                bin_names = bin_names,
            )
            
            procinfo["systematics"] = {}
            
            for systkey, systinfo in d_config["systematics"].items() :
                
                if (procname in systinfo["procs"]["bkg"]) :
                    
                    procinfo["systematics"][systkey] = systinfo
            
            if (len(procinfo["systematics"])) :
                
                get_and_set_syst_errs(
                    procname = procname,
                    inhistfile = inhistfile,
                    systs = procinfo["systematics"],
                    samples = samples,
                    scales = scales,
                    rebin = rebin,
                    bin_ids = bin_ids,
                    bin_names = bin_names,
                    hist_nom = hist_proc_nom,
                    ch_obj = cb.cp().process([procname]).era([era]),
                    target = cb,
                )
            
            procinfo["rateparams"] = {}
            
            for rpkey, rpinfo in d_config.get("rateparams", {}).items() :
                
                if (procname in rpinfo["procs"]["bkg"]) :
                    
                    procinfo["rateparams"][rpkey] = rpinfo
            
            if (len(procinfo["rateparams"])) :
                
                set_rateparams(
                    rateparams = procinfo["rateparams"],
                    ch_obj = cb.cp().process([procname]).era([era]),
                    target = cb,
                )
            
            inhistfile.Close()
    
    
    writer = ch.CardWriter(
        #"$TAG/$ANALYSIS/$MASS/channels_$CHANNEL/eras_$ERA/$ANALYSIS_$MASS_$CHANNEL_$ERA.txt",
        #"$TAG/$ANALYSIS/$MASS/channels_$CHANNEL/eras_$ERA/$ANALYSIS_$MASS_$CHANNEL_$ERA_input.root",
        "$TAG/$ANALYSIS/channels_$CHANNEL/eras_$ERA/$MASS/card_$ANALYSIS_$MASS_$CHANNEL_$ERA.txt",
        #"$TAG/$ANALYSIS/channels_$CHANNEL/eras_$ERA/$MASS/card_$ANALYSIS_$MASS_$CHANNEL_$ERA_input.root",
    )
    writer.WriteCards(args.outdir, cb)
    
    writer = ch.CardWriter(
        "$TAG/$ANALYSIS/channels_all/eras_$ERA/$MASS/card_$ANALYSIS_$MASS_$ERA.txt",
        #"$TAG/$ANALYSIS/channels_all/eras_$ERA/$MASS/card_$ANALYSIS_$MASS_$ERA_input.root",
    )
    writer.WriteCards(args.outdir, cb)
    
    writer = ch.CardWriter(
        "$TAG/$ANALYSIS/channels_all/eras_all/$MASS/card_$ANALYSIS_$MASS.txt",
        #"$TAG/$ANALYSIS/channels_all/eras_all/$MASS/card_$ANALYSIS_$MASS_input.root",
    )
    writer.WriteCards(args.outdir, cb)
    
    for ch_combo in args.chcombos :
        
        channels = ch_combo.strip().split(",")
        
        writer = ch.CardWriter(
            f"$TAG/$ANALYSIS/channels_{'_'.join(channels)}/eras_all/$MASS/card_$ANALYSIS_$MASS.txt",
            #f"$TAG/$ANALYSIS/channels_{'_'.join(channels)}/eras_all/$MASS/card_$ANALYSIS_$MASS_input.root",
        )
        
        writer.WriteCards(args.outdir, cb.cp().channel(channels))


if (__name__ == "__main__") :
    
    main()