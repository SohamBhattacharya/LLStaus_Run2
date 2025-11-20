#!/usr/bin/env python3

import argparse
import functools
import numpy
import operator
import os

from sortedcontainers import SortedList as slist

import ROOT
import CombineHarvester.CombineTools.ch as ch

import utils.commonutils as cmut
from utils.commonutils import yaml


ZERO_THRESHOLD = 0.001 # Minimum value for a bin to be considered non-zero

STAT_LABEL = "stat"
STAT_GMN_LABEL = f"{STAT_LABEL}_gmN"
STAT_GMN_LABEL_REMOVE = f"{STAT_LABEL}_gmN_REMOVE"
SYST_LABEL = "syst"

STAT_REGEX = f"{STAT_LABEL}.*"
STAT_REGEX_GMN = f"{STAT_GMN_LABEL}.*"
STAT_REGEX_NO_GMN = f"{STAT_LABEL}_(?!gmN).*"
SYST_REGEX = f"{SYST_LABEL}.*"
STAT_GMN_LABEL_REMOVE_REGEX = f"{STAT_GMN_LABEL_REMOVE}.*"

ERA_ADDED_STR = "added"


def fix_stat_gmN(l_fnames) :
    
    for fname in l_fnames :
        
        cmut.logger.info(f"Fixing gmN stat entries [{fname}]")
        nspaces = len(STAT_GMN_LABEL) - len(STAT_LABEL) - 2 # -2 for the ' 0' after gmN
        str_spaces = " "*nspaces
        #sed_cmd = f"sed -i -E \"s/{STAT_GMN_LABEL}(.+)lnN/{STAT_LABEL}\\1{str_spaces}gmN 0/g\" {fname}"
        sed_cmd = f"sed -i -E \"s/{STAT_GMN_LABEL}_([0-9]+.?[0-9]+?)_(.*)lnN/{STAT_LABEL}_\\2{str_spaces}gmN \\1/g\" {fname}"
        #print(sed_cmd)
        retval = os.system(sed_cmd)
        
        if retval :
            
            cmut.logger.error("Unable to fix.")
            exit(retval)


def set_stat_errs(
    ch_obj,
    target,
    bin_name,
    val,
    err,
    alpha = None
) :
    """
    alpha is the extrapolation factor: n=alpha*N
    n is the yield
    N is the raw count
    For e.g.: alpha can be the normalization factor xsec*lumi/Ntot
    """
    
    if val >= ZERO_THRESHOLD :
        
        err_u_rel = numpy.clip(1.0+(err/val), 0.001, 2.0)
        err_d_rel = numpy.clip(1.0-(err/val), 0.001, 2.0)
        
        #print("stat", bin_id, bin_name, val, err, err_rel)
        
        ch_obj.cp().bin([bin_name]).AddSyst(
            target = target,
            #name = "stat_$PROCESS_$BIN_chn_$CHANNEL_$ERA",
            name = f"{STAT_LABEL}_$PROCESS_$BIN",
            type = "lnN",
            valmap = ch.SystMap()(
                (err_d_rel, err_u_rel)
            )
        )
            
        # Note: ch.SystMap()((err_rel)) produces "incorrect" values when calling GetUncertainty()
        # Check description of lnN here: http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/part2/settinguptheanalysis/#a-simple-counting-experiment
        # If a single-value 1+d is provided, the relative error is evaluated as: ((1+d) - 1/(1+d))/2
        # Whereas, if 1-d/1+d is provided, the relative error is evaluated as: ((1+d) - (1-d))/2 = d
    
    # gmN error if the yield is 0
    # CombineHarvester cannot handle gmN
    # Set it as lnN first
    # Then change these lines in the card
    else :
        
        if alpha :
            
            gmN_N = int(round(val / alpha))
            
            ch_obj.cp().bin([bin_name]).AddSyst(
                target = target,
                name = f"{STAT_GMN_LABEL}_{gmN_N}_$PROCESS_$BIN",
                type = "lnN",
                valmap = ch.SystMap()(
                    (float(alpha))
                )
            )
        
        else :
            
            ch_obj.cp().bin([bin_name]).AddSyst(
                target = target,
                name = f"{STAT_LABEL}_$PROCESS_$BIN",
                type = "lnN",
                valmap = ch.SystMap()(
                    (0.001, 2.0)
                )
            )


def set_stat_errs_from_hist(
    ch_obj,
    target,
    hist,
    bin_ids,
    bin_names,
    alpha = None
) :
    """
    alpha is the extrapolation factor: n=alpha*N
    n is the yield
    N is the raw count
    For e.g.: alpha can be the normalization factor xsec*lumi/Ntot
    """
    
    for bin_id, bin_name in zip(bin_ids, bin_names) :
        
        #print(bin_id, bin_name)
        val = hist.GetBinContent(bin_id)
        val = 0 if (val < 0) else val
        
        err = hist.GetBinError(bin_id)
        
        set_stat_errs(
            ch_obj = ch_obj,
            target = target,
            bin_name = bin_name,
            val = val,
            err = err,
            alpha = alpha,
        )


def set_syst_errs(
    ch_obj,
    target,
    bin_name,
    val,
    err_u,
    err_d,
    systname,
    iseracorr,
) :
   
    # Do not allow val+-err to go too small (below 0.002)
    err_u_rel = numpy.clip(1.0+(err_u/val), 0.002/val, 2.0) if (val >= 0.002) else 1.001
    err_d_rel = numpy.clip(1.0+(err_d/val), 0.002/val, 2.0) if (val >= 0.002) else 0.999
    
    name_tmp = systname
    
    if not name_tmp.startswith(SYST_LABEL) :
        name_tmp = f"{SYST_LABEL}_{name_tmp}"
    
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

def set_syst_errs_from_hist(
    ch_obj,
    target,
    hist_nom,
    hist_u,
    hist_d,
    bin_ids,
    bin_names,
    systname,
    iseracorr,
    operation_u = None,
    operation_d = None,
) :
    
    for bin_id, bin_name in zip(bin_ids, bin_names) :
        
        val = hist_nom.GetBinContent(bin_id)
        val = 0 if (val < 0) else val
        err_u = hist_u.GetBinContent(bin_id) - val
        err_d = hist_d.GetBinContent(bin_id) - val
        
        #err_u_rel = numpy.clip(1.0+(err_u/val), 0.001, 2.0) if val else 2.0
        #err_d_rel = numpy.clip(1.0+(err_d/val), 0.001, 2.0) if val else 0.001
        
        d_fmt = {
            "u" : err_u,
            "d" : err_d,
        }
        
        err_u = eval(operation_u.format(**d_fmt)) if operation_u else err_u
        err_d = eval(operation_d.format(**d_fmt)) if operation_d else err_d
        
        #err_u_rel = numpy.clip(1.0+(err_u/val), 0.001, 2.0) if val else 1.001
        #err_d_rel = numpy.clip(1.0+(err_d/val), 0.001, 2.0) if val else 0.999
        
        ## Do not allow val+-err to go too small (below 0.002)
        #err_u_rel = numpy.clip(1.0+(err_u/val), 0.002/val, 2.0) if (val >= 0.002) else 1.001
        #err_d_rel = numpy.clip(1.0+(err_d/val), 0.002/val, 2.0) if (val >= 0.002) else 0.999
        #
        #name_tmp = f"{SYST_LABEL}_{systname}"
        #
        #if (not iseracorr) :
        #    name_tmp = f"{name_tmp}_$ERA"
        #
        #ch_obj.cp().bin([bin_name]).AddSyst(
        #    target = target,
        #    name = name_tmp,
        #    type = "lnN",
        #    valmap = ch.SystMap()(
        #        (err_d_rel, err_u_rel)
        #    )
        #)
        
        set_syst_errs(
            ch_obj = ch_obj,
            target = target,
            bin_name = bin_name,
            val = val,
            err_u = err_u,
            err_d = err_d,
            systname = systname,
            iseracorr = iseracorr
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
                
                usefile = inhistfile
                histname = systinfo[systvar]
                
                if (":" in systinfo[systvar]) :
                    
                    usefile, histname = systinfo[systvar].split(":")
                
                hist_proc_syst = cmut.get_hist(
                    histfile = usefile,
                    histname = histname,
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
                    raise Exception("Systematics not found")
                
                usefile = inhistfile
                
                if (":" in histname) :
                    
                    usefile, histname = histname.split(":")
                
                hist_proc_syst = cmut.get_hist(
                    histfile = usefile,
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
                raise Exception
            
            #hist_proc_syst.Print("range")
            d_hist_proc_syst[systname][systvar] = hist_proc_syst
        
        set_syst_errs_from_hist(
            ch_obj = ch_obj,
            target = target,
            hist_nom = hist_nom,
            hist_u = d_hist_proc_syst[systname]["u"],
            hist_d = d_hist_proc_syst[systname]["d"],
            bin_ids = bin_ids,
            bin_names = bin_names,
            systname = systname,
            iseracorr = systinfo["iseracorr"],
            operation_u = systinfo.get("operation_u", None),
            operation_d = systinfo.get("operation_d", None),
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


def eval_stat_unc(ch_obj) :
    
    cb_stat = ch_obj.cp().syst_name([STAT_REGEX_NO_GMN])
    cb_stat_gmN = ch_obj.cp().syst_name([STAT_REGEX_GMN])
    
    unc_stat_gmN = []
    cb_stat_gmN.ForEachSyst(lambda x: unc_stat_gmN.append(x.value_u()))
    
    # Use the uncertainty if yield is non-zero
    # Otherwise, add the gmN uncertainties in quadrature
    
    if ch_obj.GetRate() :
        
        unc_stat = cb_stat.GetUncertainty()
    
    else :
        
        unc_stat = -numpy.log((1-0.68)/2) * numpy.sum(numpy.array(unc_stat_gmN)**2)**0.5
    
    return float(unc_stat)


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
        help = "Custom channel combinations to create cards for; comma separated lists: ch1,ch2  ch3,ch5",
        type = str,
        nargs = "+",
        required = False,
        default = [],
    )
    
    parser.add_argument(
        "--eracombos",
        help = "Custom era combinations to create cards for; comma separated lists: era1,era2  era1,era3",
        type = str,
        nargs = "+",
        required = False,
        default = [],
    )
    
    parser.add_argument(
        "--combpars",
        help = "Create cards combining these parameters",
        type = str,
        nargs = "+",
        required = False,
        default = [],
        choices = ["channel", "era"]
    )
    
    parser.add_argument(
        "--blind",
        help = "Will not add observation to the cards; this appends \"_blinded\" to the output directory name",
        action = "store_true",
        default = False,
    )
    
    parser.add_argument(
        "--pseudodata",
        help = "Will set the observation to (total bkg), or to (signal + total bkg); this appends \"_pseudodataB\" or \"_pseudodataSB\" to the output directory name",
        type = str,
        required = False,
        default = None,
        choices = ["B", "SB"]
    )
    
    parser.add_argument(
        "--yields_uncs",
        help = "Will produce yaml files containing yields and uncertainties",
        action = "store_true",
        default = False,
    )
    
    parser.add_argument(
        "--yields_uncs_sigs",
        help = "File containing list of signals to be included in the yield and uncertainty yamls; if not provided, all signals will be included",
        type = str,
        required = False
    )
    
    parser.add_argument(
        "--sigs",
        help = "File containing list of signals to be processed; if not provided, all signals will be processed",
        type = str,
        required = False,
        default = None
    )
    
    parser.add_argument(
        "--nocards",
        help = "Will not create datacards",
        action = "store_true",
        default = False,
    )
    
    parser.add_argument(
        "--era_added_cards",
        help = "Will add (i.e. sum up eras in each bin) the eras together to create an additional set of cards",
        action = "store_true",
        default = False,
    )
    
    # Parse arguments
    args = parser.parse_args()
    #d_args = vars(args)
    
    if args.blind and not args.pseudodata:
        
        args.outdir = f"{args.outdir}_blinded"
    
    if args.pseudodata:
        
        args.outdir = f"{args.outdir}_pseudodata{args.pseudodata}"
    
    cb = ch.CombineHarvester()
    cb.SetFlag("filters-use-regex", True)
    
    l_analyses = []
    l_channels = []
    l_eras = []
    
    l_categories_era_added = []
    
    # Key structure: dict[(analysis, channel, era)]
    d_config_all = {}
    
    # Key structure: dict[analysis][channel][era]
    d_bin_ids = {}
    d_bin_names = {}
    d_sig_procs = {}
    d_bkg_procs = {}
    
    # Key structure: dict[analysis][channel][era][process]
    d_syst_names = {}
    
    # Key structure: dict[analysis][channel][era][process][systcombo]
    d_systcombo_names = {}
    
    l_bin_edges = []
    
    for cfg in args.configs :
        
        d_config = cmut.load_config(cfg)
        #print(d_config)
        
        analysis = d_config["analysis"]
        channel = d_config["channel"]
        era = d_config["era"]
        cfg_key = (analysis, channel, era)
        
        assert (cfg_key not in d_config_all), f"(analysis, channel, era) = {cfg_key} is repeated. Check input configurations."
        d_config_all[cfg_key] = d_config
        
        if analysis not in d_bin_ids :
            
            d_bin_ids[analysis] = {}
            d_bin_names[analysis] = {}
            d_sig_procs[analysis] = {}
            d_bkg_procs[analysis] = {}
            d_syst_names[analysis] = {}
            d_systcombo_names[analysis] = {}
        
        if channel not in d_bin_ids[analysis] :
            
            d_bin_ids[analysis][channel] = {}
            d_bin_names[analysis][channel] = {}
            d_sig_procs[analysis][channel] = {}
            d_bkg_procs[analysis][channel] = {}
            d_syst_names[analysis][channel] = {}
            d_systcombo_names[analysis][channel] = {}
        
        if era not in d_bin_ids[analysis][channel] :
            
            d_bin_ids[analysis][channel][era] = []
            d_bin_names[analysis][channel][era] = []
            d_sig_procs[analysis][channel][era] = []
            d_bkg_procs[analysis][channel][era] = []
            d_syst_names[analysis][channel][era] = {}
            d_systcombo_names[analysis][channel][era] = {}
        
        l_analyses.append(analysis)
        l_channels.append(channel)
        l_eras.append(era)
        
        # Copy the configs to the output directory
        outdir = f"{args.outdir}/{analysis}"
        outdir_configs = f"{outdir}/configs"
        os.system(f"mkdir -p {outdir}")
        os.system(f"mkdir -p {outdir_configs}")
        cfg_name, cfg_ext = os.path.splitext(cfg)
        cfg_cpy = f"{outdir_configs}/{os.path.basename(cfg_name)}_ana_{analysis}_chn_{channel}_era_{era}{cfg_ext}"
        os.system(f"cp -vf {cfg} {cfg_cpy}")
    
    
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
        channel = d_config["channel"]
        era = d_config["era"]
        lumi = eval(d_config["lumi"]) if isinstance(d_config["lumi"], str) else d_config["lumi"]
        binstart = d_config["binstart"]
        categories = [(_binnum, f"bin{_binnum}_{channel}_{era}") for _binnum in range(binstart, binstart+nbins)]
        l_categories_era_added.extend([(_binnum, f"bin{_binnum}_{channel}_{ERA_ADDED_STR}") for _binnum in range(binstart, binstart+nbins)])
        bin_ids = [_cat[0] for _cat in categories]
        bin_names = [_cat[1] for _cat in categories]
        l_sig_processes = []
        l_sig_processes_yields_uncs = []
        l_bkg_processes = list(d_config["bkg"]["procs"].keys())
        
        d_bin_ids[analysis][channel][era].extend(bin_ids)
        d_bin_names[analysis][channel][era].extend(bin_names)
        
        print(bin_ids, bin_names)
        
        if ("sig" in d_config) :
            
            if ("loadconfig" in d_config["sig"]) :
                
                d_config["sig"].update(cmut.load_config(d_config["sig"]["loadconfig"]))
            
            l_sig_processes = list(d_config["sig"]["procs"].keys())
            
            if args.sigs :
                l_sig_processes_tmp = numpy.loadtxt(args.sigs, dtype = str)
                l_sig_processes = [_proc for _proc in l_sig_processes if _proc in l_sig_processes_tmp]
            
            if args.yields_uncs_sigs :
                l_sig_processes_yields_uncs = numpy.loadtxt(args.yields_uncs_sigs, dtype = str)
            
            # If --nocards is set, then only include the signal processes for which the yield and uncertainty yamls are to be stored
            if args.nocards and args.yields_uncs_sigs :
                
                l_sig_processes = [_proc for _proc in l_sig_processes if _proc in l_sig_processes_yields_uncs]
            
            cb.AddProcesses(
                mass = l_sig_processes,
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
        
        d_xsec = cmut.load_config(d_config["bkg"]["xsecfile"])
        #print(d_xsec)
        
        if ("sig" in d_config) :
            
            inhistfile_sig = ROOT.TFile.Open(d_config["sig"]["histfile"])
            cutflows_sig = cmut.load_config(d_config["sig"]["cutflowsfile"])
            
            scale_allprocs_sig = eval(d_config["sig"].get("scaleby_allprocs", "1.0"))
            
            for procname, procinfo in d_config["sig"]["procs"].items() :
                
                if procname not in l_sig_processes :
                    
                    continue
                
                samples = procinfo["samples"]
                issusy = procinfo["issusy"]
                neventkey = d_config["sig"]["neventkey"]
                neventkey_lifetime_reweighting = d_config["sig"]["neventkey_lifetime_reweighting"] if issusy else None
                ismc = True if issusy else procinfo["ismc"]
                
                d_sig_procs[analysis][channel][era].append(procname)
                
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
                l_nevents_sample = []
                
                for sample in samples :
                    
                    scale = eval(procinfo["scaleby"]) * scale_allprocs_sig
                    
                    if (ismc) :
                        
                        scale *= lumi
                        
                        if (procinfo["xsnorm"]) :
                            
                            if (issusy and "_from_ctau" in sample) :
                                
                                neventtot = functools.reduce(operator.getitem, [sample]+neventkey_lifetime_reweighting.split("."), cutflows_sig)
                            
                            else :
                                
                                neventtot = functools.reduce(operator.getitem, [sample]+neventkey.split("."), cutflows_sig)
                            
                            scale *= (d_xsec[sample] / neventtot)
                            
                            l_nevents_sample.append(neventtot)
                    
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
                
                if not l_bin_edges :
                    
                    l_bin_edges = [hist_proc_nom.GetXaxis().GetBinLowEdge(_bin) for _bin in bin_ids] + [hist_proc_nom.GetXaxis().GetBinUpEdge(bin_ids[-1])]
                
                #hist_proc_nom.Print("range")
                
                alpha = procinfo.get(
                    "alpha",
                    numpy.sum(numpy.power(scales, 2))**0.5
                )
                
                cb.cp().mass([procname]).channel([channel]).process(["sig"]).era([era]).ForEachProc(
                    lambda x : x.set_rate(
                        #hist_proc_nom.GetBinContent(x.bin_id()) if hist_proc_nom.GetBinContent(x.bin_id()) else 0.001
                        #hist_proc_nom.GetBinContent(x.bin_id())
                        0.001 if (not hist_proc_nom.GetBinContent(x.bin_id()) and not alpha) else hist_proc_nom.GetBinContent(x.bin_id())
                ))
                
                set_stat_errs_from_hist(
                    ch_obj = cb.cp().mass([procname]).process(["sig"]).era([era]),
                    target = cb,
                    hist = hist_proc_nom,
                    bin_ids = bin_ids,
                    bin_names = bin_names,
                    alpha = alpha,
                )
                
                procinfo["systematics"] = {}
                
                for systkey, systinfo in d_config["systematics"].items() :
                    
                    if (procname in systinfo["procs"]["sig"] or "sig" in systinfo["procs"]["sig"]) :
                        
                        procinfo["systematics"][systkey] = systinfo
                
                d_systcombo_names[analysis][channel][era][procname] = {}
                
                for systkey, systinfo in d_config.get("systcombos", {}).items() :
                    
                    if (procname in systinfo["procs"]["sig"] or "sig" in systinfo["procs"]["sig"]) :
                        
                        d_systcombo_names[analysis][channel][era][procname][systkey] = systinfo["names"]
                
                if (len(procinfo["systematics"])) :
                    
                    #l_systematics.extend(procinfo["systematics"].keys())
                    d_syst_names[analysis][channel][era][procname] = list(procinfo["systematics"].keys())
                    
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
        
        
        scale_allprocs_bkg = eval(d_config["bkg"].get("scaleby_allprocs", "1.0"))
        
        for procname, procinfo in d_config["bkg"]["procs"].items() :
            
            ismc = procinfo["ismc"]
            inhistfile = ROOT.TFile.Open(procinfo["histfile"])
            cutflows_bkg = cmut.load_config(d_config["bkg"]["cutflowsfile"])
            neventkey = d_config["bkg"]["neventkey"]
            samples = procinfo["samples"]
            
            d_bkg_procs[analysis][channel][era].append(procname)
            
            scales = []
            l_nevents_sample = []
            
            for sample in samples :
                
                scale = eval(procinfo["scaleby"]) * scale_allprocs_bkg
                
                if (ismc) :
                    
                    scale *= lumi
                    
                    if (procinfo["xsnorm"]) :
                        
                        neventtot = functools.reduce(operator.getitem, [sample]+neventkey.split("."), cutflows_bkg)
                        scale *= (d_xsec[sample] / neventtot)
                        
                        l_nevents_sample.append(neventtot)
                
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
            
            if not l_bin_edges :
                
                l_bin_edges = [hist_proc_nom.GetXaxis().GetBinLowEdge(_bin) for _bin in bin_ids] + [hist_proc_nom.GetXaxis().GetBinUpEdge(bin_ids[-1])]
            
            alpha = procinfo.get(
                "alpha",
                numpy.sum(numpy.power(scales, 2))**0.5
            )
            
            cb.cp().channel([channel]).process([procname]).era([era]).ForEachProc(
                lambda x : x.set_rate(
                    #hist_proc_nom.GetBinContent(x.bin_id()) if hist_proc_nom.GetBinContent(x.bin_id()) else 0.001
                    #hist_proc_nom.GetBinContent(x.bin_id())
                    0.001 if (not hist_proc_nom.GetBinContent(x.bin_id()) and not alpha) else hist_proc_nom.GetBinContent(x.bin_id())
            ))
            
            print(f"[{procname}] [alpha {alpha}]")
            set_stat_errs_from_hist(
                ch_obj = cb.cp().process([procname]).era([era]),
                target = cb,
                hist = hist_proc_nom,
                bin_ids = bin_ids,
                bin_names = bin_names,
                alpha = alpha
            )
            
            procinfo["systematics"] = {}
            
            for systkey, systinfo in d_config["systematics"].items() :
                
                if (procname in systinfo["procs"]["bkg"]) :
                    
                    procinfo["systematics"][systkey] = systinfo
            
            d_systcombo_names[analysis][channel][era][procname] = {}
            
            for systkey, systinfo in d_config.get("systcombos", {}).items() :
                
                if (procname in systinfo["procs"]["bkg"]) :
                    
                    d_systcombo_names[analysis][channel][era][procname][systkey] = systinfo["names"]
            
            if (len(procinfo["systematics"])) :
                
                #l_systematics.extend(procinfo["systematics"].keys())
                d_syst_names[analysis][channel][era][procname] = list(procinfo["systematics"].keys())
                
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
        
        load_obs = ((not args.blind) and ("obs" in d_config) and (not args.pseudodata))
        
        if (load_obs or args.pseudodata) :
            
            if (load_obs or args.pseudodata == "B") :
                
                cb.AddObservations(
                    mass = ["*"],
                    analysis = [analysis],
                    era = [era],
                    channel = [channel],
                    bin = categories,
                )
            
            elif (args.pseudodata == "SB") :
                
                cb.AddObservations(
                    mass = l_sig_processes,
                    analysis = [analysis],
                    era = [era],
                    channel = [channel],
                    bin = categories,
                )
            
            if (not args.pseudodata) :
                
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
                
                inhistfile_obs.Close()
        
        if (args.pseudodata == "B") :
            
            cb.cp().channel([channel]).era([era]).ForEachObs(
                lambda x : x.set_rate(
                    cb.cp().channel([channel]).process(l_bkg_processes).era([era]).bin_id([x.bin_id()]).GetRate()
            ))
        
        elif (args.pseudodata == "SB") :
            
            for procname in l_sig_processes :
                
                cb.cp().mass([procname]).channel([channel]).era([era]).ForEachObs(
                    lambda x : x.set_rate(
                        cb.cp().channel([channel]).process(l_bkg_processes).era([era]).bin_id([x.bin_id()]).GetRate() +
                        cb.cp().mass([procname]).channel([channel]).process(["sig"]).era([era]).bin_id([x.bin_id()]).GetRate()
                ))
    
    if args.yields_uncs :
        
        for analysis in slist(d_bin_ids.keys()) :
            
            for channel in slist(d_bin_ids[analysis].keys()) :
                
                # Each era
                l_eras_tmp = [[_era] for _era in d_bin_ids[analysis][channel].keys()]
                
                # Also the combination of all eras
                if len(d_bin_ids[analysis][channel].keys()) > 1 :
                    
                    l_eras_tmp.append(list(d_bin_ids[analysis][channel].keys()))
                
                for ieras, eras in enumerate(l_eras_tmp) :
                    
                    d_yields = {}
                    d_systematics = {}
                    
                    eras_str = ".".join(eras)
                    d_yields[eras_str] = {}
                    d_systematics[eras_str] = {}
                    
                    l_bin_ids = []
                    l_sig_processes = []
                    l_bkg_processes = []
                    
                    for era in eras :
                        
                        l_bin_ids.extend(d_bin_ids[analysis][channel][era])
                        l_sig_processes.extend(d_sig_procs[analysis][channel][era])
                        l_bkg_processes.extend(d_bkg_procs[analysis][channel][era])
                    
                    # Remove duplicate entries
                    l_bin_ids = slist(set(l_bin_ids))
                    l_sig_processes = slist(set(l_sig_processes))
                    l_bkg_processes = slist(set(l_bkg_processes))
                    
                    if args.yields_uncs_sigs :
                        
                        #l_sig_processes_sel = numpy.loadtxt(args.yields_uncs_sigs, dtype = str)
                        #l_sig_processes = [_proc for _proc in l_sig_processes if _proc in l_sig_processes_sel]
                        l_sig_processes = [_proc for _proc in l_sig_processes if _proc in l_sig_processes_yields_uncs]
                    
                    if not args.blind :
                        
                        d_yields[eras_str]["obs"] = {}
                        
                        for bin_id in l_bin_ids :
                            
                            d_yields[eras_str]["obs"][bin_id] = {
                                "yield": cb.cp().era(eras).bin_id([bin_id]).GetObservedRate()
                            }
                    
                    for bkg_proc in l_bkg_processes :
                        
                        l_systematics = []
                        d_systcombos = {}
                        
                        for era in eras :
                            
                            l_systematics.extend(d_syst_names[analysis][channel][era][bkg_proc])
                            
                            for syst_name, combolist in d_systcombo_names[analysis][channel][era][bkg_proc].items() :
                                
                                if (syst_name not in d_systcombos) :
                                    d_systcombos[syst_name] = []
                                
                                d_systcombos[syst_name].extend(combolist)
                        
                        # Remove duplicates
                        l_systematics = slist(set(l_systematics))
                        
                        for syst_name, combolist in d_systcombos.items() :
                            
                            d_systcombos[syst_name] = set(combolist)
                        
                        # Per bkg, per bin yields and systematics
                        d_yields[eras_str][bkg_proc] = {}
                        d_systematics[eras_str][bkg_proc] = {_syst: {} for _syst in l_systematics+list(d_systcombos.keys())}
                        
                        print(eras, bkg_proc, l_bin_ids, l_systematics)
                        
                        for bin_id in l_bin_ids :
                            
                            cb_bkg = cb.cp().era(eras).process([bkg_proc]).bin_id([bin_id])
                            #cb_bkg_stat = cb.cp().era(eras).process([bkg_proc]).bin_id([bin_id]).syst_name([STAT_REGEX_NO_GMN])
                            #cb_bkg_stat_gmN = cb.cp().era(eras).process([bkg_proc]).bin_id([bin_id]).syst_name([STAT_REGEX_GMN])
                            cb_bkg_syst = cb.cp().era(eras).process([bkg_proc]).bin_id([bin_id]).syst_name([SYST_REGEX])
                            
                            #unc_stat_gmN = []
                            #cb_bkg_stat_gmN.ForEachSyst(lambda x: unc_stat_gmN.append(x.value_u()))
                            #
                            # Use the uncertainty if yield is non-zero
                            # Otherwise, add the gmN uncertainties in quadrature
                            #unc_stat = cb_bkg_stat.GetUncertainty() if cb_bkg.GetRate() else -numpy.log((1-0.68)/2)*sum(numpy.array(unc_stat_gmN)**2)**0.5
                            
                            d_yields[eras_str][bkg_proc][bin_id] = {
                                "yield": cb_bkg.GetRate(),
                                "unc_stat": eval_stat_unc(ch_obj = cb_bkg.cp()),
                                "unc_syst": cb_bkg_syst.GetUncertainty(),
                            }
                            
                            #d_systematics[eras_str][bkg_proc][bin_id] = {}
                            
                            for syst_name in l_systematics :
                                
                                syst_rgx = f"{SYST_LABEL}_{syst_name}.*"
                                #print(syst_rgx)
                                cb_bkg_syst = cb.cp().era(eras).process([bkg_proc]).bin_id([bin_id]).syst_name([syst_rgx])
                                
                                d_systematics[eras_str][bkg_proc][syst_name][bin_id] = {
                                    "yield": cb_bkg.GetRate(),
                                    "unc_abs": cb_bkg_syst.GetUncertainty(),
                                    "unc_rel": cb_bkg_syst.GetUncertainty()/cb_bkg.GetRate() if cb_bkg.GetRate() else -1.0,
                                }
                            
                            for syst_name, combolist in d_systcombos.items() :
                                
                                l_syst_rgx = [f"{SYST_LABEL}_{_syst}.*" for _syst in combolist]
                                cb_bkg_syst = cb.cp().era(eras).process([bkg_proc]).bin_id([bin_id]).syst_name(l_syst_rgx)
                                
                                d_systematics[eras_str][bkg_proc][syst_name][bin_id] = {
                                    "yield": cb_bkg.GetRate(),
                                    "unc_abs": cb_bkg_syst.GetUncertainty(),
                                    "unc_rel": cb_bkg_syst.GetUncertainty()/cb_bkg.GetRate() if cb_bkg.GetRate() else -1.0,
                                }
                        
                        # Per bkg, all bin yields
                        cb_bkg = cb.cp().era(eras).process([bkg_proc])
                        #cb_bkg_stat = cb.cp().era(eras).process([bkg_proc]).syst_name([STAT_REGEX_NO_GMN])
                        cb_bkg_syst = cb.cp().era(eras).process([bkg_proc]).syst_name([SYST_REGEX])
                        
                        d_yields[eras_str][bkg_proc]["total"] = {
                            "yield": cb_bkg.GetRate(),
                            #"unc_stat": cb_bkg_stat.GetUncertainty(),
                            "unc_stat": eval_stat_unc(ch_obj = cb_bkg.cp()),
                            "unc_syst": cb_bkg_syst.GetUncertainty(),
                        }
                    
                    if not args.blind :
                        
                        d_yields[eras_str]["obs"]["total"] = {
                            "yield": cb.cp().era(eras).GetObservedRate()
                        }
                    
                    # All bkg, per bin yields
                    d_yields[eras_str]["total_bkg"] = {}
                    
                    for bin_id in l_bin_ids :
                        
                        cb_bkg = cb.cp().era(eras).process(l_bkg_processes).bin_id([bin_id])
                        #cb_bkg_stat = cb.cp().era(eras).process(l_bkg_processes).bin_id([bin_id]).syst_name([STAT_REGEX_NO_GMN])
                        cb_bkg_syst = cb.cp().era(eras).process(l_bkg_processes).bin_id([bin_id]).syst_name([SYST_REGEX])
                        
                        d_yields[eras_str]["total_bkg"][bin_id] = {
                            "yield": cb_bkg.GetRate(),
                            #"unc_stat": cb_bkg_stat.GetUncertainty(),
                            "unc_stat": eval_stat_unc(ch_obj = cb_bkg.cp()),
                            "unc_syst": cb_bkg_syst.GetUncertainty(),
                        }
                    
                    # All bkg, all bin yields
                    cb_bkg = cb.cp().era(eras).process(l_bkg_processes)
                    #cb_bkg_stat = cb.cp().era(eras).process(l_bkg_processes).syst_name([STAT_REGEX_NO_GMN])
                    cb_bkg_syst = cb.cp().era(eras).process(l_bkg_processes).syst_name([SYST_REGEX])
                    
                    d_yields[eras_str]["total_bkg"]["total"] = {
                        "yield": cb_bkg.GetRate(),
                        #"unc_stat": cb_bkg_stat.GetUncertainty(),
                        "unc_stat": eval_stat_unc(ch_obj = cb_bkg.cp()),
                        "unc_syst": cb_bkg_syst.GetUncertainty(),
                    }
                    
                    for sig_proc in l_sig_processes :
                        
                        l_systematics = []
                        d_systcombos = {}
                        
                        for era in eras :
                            
                            l_systematics.extend(d_syst_names[analysis][channel][era][sig_proc])
                            
                            for syst_name, combolist in d_systcombo_names[analysis][channel][era][sig_proc].items() :
                                
                                if (syst_name not in d_systcombos) :
                                    d_systcombos[syst_name] = []
                                
                                d_systcombos[syst_name].extend(combolist)
                        
                        # Remove duplicates
                        l_systematics = slist(set(l_systematics))
                        
                        for syst_name, combolist in d_systcombos.items() :
                            
                            d_systcombos[syst_name] = set(combolist)
                        
                        # Per sig, per bin yields and systematics
                        d_yields[eras_str][sig_proc] = {}
                        d_systematics[eras_str][sig_proc] = {}
                        d_systematics[eras_str][sig_proc] = {_syst: {} for _syst in l_systematics+list(d_systcombos.keys())}
                        
                        print(eras, sig_proc, l_bin_ids, l_systematics)
                        
                        for bin_id in l_bin_ids :
                            
                            #print(sig_proc, bin_id)
                            cb_sig = cb.cp().era(eras).mass([sig_proc]).bin_id([bin_id])
                            #cb_sig_stat = cb.cp().era(eras).mass([sig_proc]).bin_id([bin_id]).syst_name([STAT_REGEX_NO_GMN])
                            cb_sig_syst = cb.cp().era(eras).mass([sig_proc]).bin_id([bin_id]).syst_name([SYST_REGEX])
                            
                            d_yields[eras_str][sig_proc][bin_id] = {
                                "yield": cb_sig.GetRate(),
                                #"unc_stat": cb_sig_stat.GetUncertainty(),
                                "unc_stat": eval_stat_unc(ch_obj = cb_sig.cp()),
                                "unc_syst": cb_sig_syst.GetUncertainty(),
                            }
                            
                            #d_systematics[eras_str][sig_proc][bin_id] = {}
                            
                            for syst_name in l_systematics :
                                
                                syst_rgx = f"{SYST_LABEL}_{syst_name}.*"
                                cb_sig_syst = cb.cp().era(eras).mass([sig_proc]).bin_id([bin_id]).syst_name([syst_rgx])
                                
                                d_systematics[eras_str][sig_proc][syst_name][bin_id] = {
                                    "yield": cb_sig.GetRate(),
                                    "unc_abs": cb_sig_syst.GetUncertainty(),
                                    "unc_rel": cb_sig_syst.GetUncertainty()/cb_sig.GetRate() if cb_sig.GetRate() else -1.0,
                                }
                            
                            for syst_name, combolist in d_systcombos.items() :
                                
                                l_syst_rgx = [f"{SYST_LABEL}_{_syst}.*" for _syst in combolist]
                                #print(syst_name, l_syst_rgx)
                                cb_sig_syst = cb.cp().era(eras).mass([sig_proc]).bin_id([bin_id]).syst_name(l_syst_rgx)
                                
                                d_systematics[eras_str][sig_proc][syst_name][bin_id] = {
                                    "yield": cb_sig.GetRate(),
                                    "unc_abs": cb_sig_syst.GetUncertainty(),
                                    "unc_rel": cb_sig_syst.GetUncertainty()/cb_sig.GetRate() if cb_sig.GetRate() else -1.0,
                                }
                        
                        # Per sig, all bin yields
                        cb_sig = cb.cp().era(eras).mass([sig_proc])
                        #cb_sig_stat = cb.cp().era(eras).mass([sig_proc]).syst_name([STAT_REGEX_NO_GMN])
                        cb_sig_syst = cb.cp().era(eras).mass([sig_proc]).syst_name([SYST_REGEX])
                        
                        d_yields[eras_str][sig_proc]["total"] = {
                            "yield": cb_sig.GetRate(),
                            #"unc_stat": cb_sig_stat.GetUncertainty(),
                            "unc_stat": eval_stat_unc(ch_obj = cb_sig.cp()),
                            "unc_syst": cb_sig_syst.GetUncertainty(),
                        }
                    
                    d_yields[eras_str]["bin_edges"] = l_bin_edges
                    d_yields[eras_str]["binstart"] = binstart
                    d_yields[eras_str]["nbins"] = nbins
                    
                    outdir = f"{args.outdir}/{analysis}/yields_and_systematics"
                    os.system(f"mkdir -p {outdir}")
                    
                    outyamlname = f"{outdir}/yields_channels_{channel}_eras_{eras_str}.yaml"
                    with open(outyamlname, "w", encoding = "utf-8") as fopen:
                        
                        yaml.dump(d_yields, fopen)
                    
                    cmut.logger.info(f"Written yields [{outyamlname}]")
                    
                    outyamlname = f"{outdir}/systematics_channels_{channel}_eras_{eras_str}.yaml"
                    with open(outyamlname, "w", encoding = "utf-8") as fopen:
                        
                        yaml.dump(d_systematics, fopen)
                    
                    cmut.logger.info(f"Written systematics [{outyamlname}]")
                    
                    # Era added card
                    if args.era_added_cards and len(l_eras) > 1 :
                        
                        l_categories_era_added = slist(set(l_categories_era_added))
                        
                        if ieras == len(l_eras_tmp)-1 :
                            
                            cb_era_added = ch.CombineHarvester()
                            cb_era_added.SetFlag("filters-use-regex", True)
                            
                            cb_era_added.AddProcesses(
                                mass = l_sig_processes,
                                analysis = [analysis],
                                era = [ERA_ADDED_STR],
                                channel = [channel],
                                procs = ["sig"],
                                bin = l_categories_era_added,
                                signal = True,
                            )
                            
                            cb_era_added.AddProcesses(
                                mass = ["*"],
                                analysis = [analysis],
                                era = [ERA_ADDED_STR],
                                channel = [channel],
                                procs = l_bkg_processes,
                                bin = l_categories_era_added,
                                signal = False,
                            )
                            
                            cb_era_added.AddObservations(
                                mass = ["*"],
                                analysis = [analysis],
                                era = [ERA_ADDED_STR],
                                channel = [channel],
                                bin = l_categories_era_added,
                            )
                            
                            for bkg_proc in l_bkg_processes :
                                
                                cb_tmp = cb.cp().channel([channel]).process([bkg_proc])
                                cb_era_added_tmp = cb_era_added.cp().era([ERA_ADDED_STR]).channel([channel]).process([bkg_proc])
                                
                                cb_era_added_tmp.ForEachProc(
                                    lambda x : x.set_rate(
                                        d_yields[eras_str][bkg_proc][x.bin_id()]["yield"]
                                ))
                                
                                l_systematics = cb_tmp.cp().syst_name([SYST_REGEX]).syst_name_set()
                                #print("*"*10, l_systematics)
                                
                                for bin_id, bin_name in l_categories_era_added :
                                    
                                    set_stat_errs(
                                        ch_obj = cb_era_added_tmp,
                                        target = cb_era_added,
                                        bin_name = bin_name,
                                        val = d_yields[eras_str][bkg_proc][bin_id]["yield"],
                                        err = d_yields[eras_str][bkg_proc][bin_id]["unc_stat"],
                                        alpha = None
                                    )
                                    
                                    for systname in l_systematics :
                                        
                                        cb_tmp_syst = cb_tmp.cp().bin_id([bin_id]).syst_name([systname])
                                        
                                        set_syst_errs(
                                            ch_obj = cb_era_added_tmp,
                                            target = cb_era_added,
                                            bin_name = bin_name,
                                            val = d_yields[eras_str][bkg_proc][bin_id]["yield"],
                                            err_u = cb_tmp_syst.GetUncertainty(),
                                            err_d = -cb_tmp_syst.GetUncertainty(),
                                            systname = systname,
                                            iseracorr = True,
                                        )
                            
                            #print("$"*20, l_sig_processes)
                            for sig_proc in l_sig_processes :
                                
                                cb_tmp = cb.cp().channel([channel]).mass([sig_proc])
                                cb_era_added_tmp = cb_era_added.cp().era([ERA_ADDED_STR]).channel([channel]).mass([sig_proc])
                                
                                cb_era_added_tmp.ForEachProc(
                                    lambda x : x.set_rate(
                                        d_yields[eras_str][sig_proc][x.bin_id()]["yield"]
                                ))
                                
                                l_systematics = cb_tmp.cp().syst_name([SYST_REGEX]).syst_name_set()
                                
                                for bin_id, bin_name in l_categories_era_added :
                                    
                                    set_stat_errs(
                                        ch_obj = cb_era_added_tmp,
                                        target = cb_era_added,
                                        bin_name = bin_name,
                                        val = d_yields[eras_str][sig_proc][bin_id]["yield"],
                                        err = d_yields[eras_str][sig_proc][bin_id]["unc_stat"],
                                        alpha = None
                                    )
                                    
                                    for systname in l_systematics :
                                        
                                        cb_tmp_syst = cb_tmp.cp().bin_id([bin_id]).syst_name([systname])
                                        
                                        set_syst_errs(
                                            ch_obj = cb_era_added_tmp,
                                            target = cb_era_added,
                                            bin_name = bin_name,
                                            val = d_yields[eras_str][sig_proc][bin_id]["yield"],
                                            err_u = cb_tmp_syst.GetUncertainty(),
                                            err_d = -cb_tmp_syst.GetUncertainty(),
                                            systname = systname,
                                            iseracorr = False,
                                        )
                            
                            cb_era_added.cp().era([ERA_ADDED_STR]).channel([channel]).ForEachObs(
                                lambda x : x.set_rate(
                                    d_yields[eras_str]["obs"][x.bin_id()]["yield"]
                            ))
                            
                            cb_era_added.SetFlag("filters-use-regex", False)
                            #cb_era_added.PrintAll()
                            
                            writer = ch.CardWriter(
                                "$TAG/$ANALYSIS/channels_$CHANNEL/eras_added/$MASS/card_ana_$ANALYSIS_mss_$MASS_chn_$CHANNEL_era_$ERA.txt",
                                "/tmp/dummy_combineharvester_input.root",
                            )
                            writer_ret = writer.WriteCards(args.outdir, cb_era_added)
                            fix_stat_gmN([_key for _key, _val in writer_ret])
    
    
    if not args.nocards :
        
        # The CardWriter complains otherwise
        cb.SetFlag("filters-use-regex", False)
        
        # *_input.root does not contain anything
        # Set it to a dummy file
        
        # Per analysis, channel, era, mass
        writer = ch.CardWriter(
            #"$TAG/$ANALYSIS/$MASS/channels_$CHANNEL/eras_$ERA/$ANALYSIS_mass_$MASS_chn_$CHANNEL_$ERA.txt",
            #"$TAG/$ANALYSIS/$MASS/channels_$CHANNEL/eras_$ERA/$ANALYSIS_mass_$MASS_chn_$CHANNEL_$ERA_input.root",
            "$TAG/$ANALYSIS/channels_$CHANNEL/eras_$ERA/$MASS/card_ana_$ANALYSIS_mss_$MASS_chn_$CHANNEL_era_$ERA.txt",
            #"$TAG/$ANALYSIS/channels_$CHANNEL/eras_$ERA/$MASS/card_ana_$ANALYSIS_mss_$MASS_chn_$CHANNEL_$ERA_era_input.root",
            "/tmp/dummy_combineharvester_input.root",
        )
        writer_ret = writer.WriteCards(args.outdir, cb)
        fix_stat_gmN([_key for _key, _val in writer_ret])
        
        if ("channel" in args.combpars) :
            
            # Per analysis, era, mass
            writer = ch.CardWriter(
                "$TAG/$ANALYSIS/channels_all/eras_$ERA/$MASS/card_ana_$ANALYSIS_mss_$MASS_chn_all_era_$ERA.txt",
                #"$TAG/$ANALYSIS/channels_all/eras_$ERA/$MASS/card_ana_$ANALYSIS_mss_$MASS_chn_all_$ERA_era_input.root",
                "/tmp/dummy_combineharvester_input.root",
            )
            writer_ret = writer.WriteCards(args.outdir, cb)
            fix_stat_gmN([_key for _key, _val in writer_ret])
        
        if ("era" in args.combpars) :
            
            # Per analysis, channel, mass
            writer = ch.CardWriter(
                "$TAG/$ANALYSIS/channels_$CHANNEL/eras_all/$MASS/card_ana_$ANALYSIS_mss_$MASS_chn_$CHANNEL_era_all.txt",
                #"$TAG/$ANALYSIS/channels_$CHANNEL/eras_all/$MASS/card_ana_$ANALYSIS_mss_$MASS_chn_$CHANNEL_era_all_input.root",
                "/tmp/dummy_combineharvester_input.root",
            )
            writer_ret = writer.WriteCards(args.outdir, cb)
            fix_stat_gmN([_key for _key, _val in writer_ret])
        
        if ("channel" in args.combpars and "era" in args.combpars) :
            
            # Per analysis, mass
            writer = ch.CardWriter(
                "$TAG/$ANALYSIS/channels_all/eras_all/$MASS/card_ana_$ANALYSIS_mss_$MASS_chn_all_era_all.txt",
                #"$TAG/$ANALYSIS/channels_all/eras_all/$MASS/card_ana_$ANALYSIS_mss_$MASS_chn_all_era_all_input.root",
                "/tmp/dummy_combineharvester_input.root",
            )
            writer_ret = writer.WriteCards(args.outdir, cb)
            fix_stat_gmN([_key for _key, _val in writer_ret])
        
        for ch_combo in args.chcombos :
            
            channels = ch_combo.strip().split(",")
            channels_tag = ".".join(channels)
            
            if not set(channels).issubset(l_channels) :
                
                cmut.logger.error(f"Invalid era in --chcombos: {args.chcombos}")
                raise Exception
            
            writer = ch.CardWriter(
                f"$TAG/$ANALYSIS/channels_{channels_tag}/eras_all/$MASS/card_ana_$ANALYSIS_mss_$MASS_chn_{channels_tag}_era_all.txt",
                "/tmp/dummy_combineharvester_input.root",
            )
            
            writer_ret = writer.WriteCards(args.outdir, cb.cp().channel(channels))
            fix_stat_gmN([_key for _key, _val in writer_ret])
        
        for era_combo in args.eracombos :
            
            eras = era_combo.strip().split(",")
            eras_tag = ".".join(eras)
            
            if not set(eras).issubset(l_eras) :
                
                cmut.logger.error(f"Invalid era in --eracombos: {args.eracombos}")
                raise Exception
            
            writer = ch.CardWriter(
                f"$TAG/$ANALYSIS/channels_all/eras_{eras_tag}/$MASS/card_ana_$ANALYSIS_mss_$MASS_chn_all_era_{eras_tag}.txt",
                "/tmp/dummy_combineharvester_input.root",
            )
            
            writer_ret = writer.WriteCards(args.outdir, cb.cp().era(eras))
            fix_stat_gmN([_key for _key, _val in writer_ret])


if (__name__ == "__main__") :
    
    main()