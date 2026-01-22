#!/usr/bin/env python3

import argparse
import itertools
import functools
from matplotlib.pylab import f
from matplotlib.pyplot import sca
import numpy
import operator
import os

import utils.commonutils as cmut
from utils.commonutils import yaml

def main() :
    
    # Argument parser
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        "--config",
        help = "Configuration yaml file.",
        type = str,
        required = True,
    )
    
    #parser.add_argument(
    #    "--tex",
    #    help = "Latex template file. Will use this to create a pdf and of the table. Pass empty string \"\" to disable latex compilation",
    #    type = str,
    #    required = False,
    #    default = "utils/tex/template_standalone_png_no-CMS.tex",
    #)
    
    parser.add_argument(
        "--outdir",
        help = "Output directory",
        type = str,
        required = False,
        default = "tmp/test_cutflows"
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    d_config = cmut.load_config(args.config)
    analysis = d_config["analysis"]
    channel = d_config["channel"]
    era_str = ".".join([_erainfo["label"] for _erainfo in d_config["eras"]])
    
    outfname = f"{args.outdir}/{analysis}/{channel}/{era_str}/cutflows_{analysis}_{channel}_{era_str}.yaml"
    outdir = os.path.dirname(outfname)
    os.system(f"mkdir -p {outdir}")
    
    if isinstance(d_config["procs"], str) :
        d_config["procs"] = cmut.load_config(d_config["procs"])["procs"]
    
    #print(d_config)
    
    d_xsec = {}
    
    d_xsec.update(cmut.get_stau_xsec_dict(
        l_samplestr = list(itertools.chain.from_iterable([_proc["samples"] for _proc in d_config["procs"].values()])),
        xsecfile = d_config["xsecfile"],
        regexp = "SMS-TStauStau_MStau-(?P<mstau>\d+)_ctau-(?P<ctau>\w+)_mLSP-(?P<mlsp>\d+)",
    ))
    #print(d_xsec)
    
    d_cutflow_input = {}
    for d_erainfo in d_config["eras"] :
        
        d_cutflow_input[d_erainfo["label"]] = cmut.load_config(d_erainfo["cutflowsfile"])
    
    d_cutflow_output = {}
    
    for prockey, d_procinfo in d_config["procs"].items() :
        
        cmut.logger.info(f"Preparing cutflow [proc {prockey}]")
        
        proclabel = d_procinfo.get("label", prockey)
        
        if d_procinfo["issusy"] :
            
            d_parse_result = cmut.parse_string_regex(
                prockey,
                regexp = "SMS-TStauStau_MStau-(?P<mstau>\d+)_ctau-(?P<ctau>\w+)_mLSP-(?P<mlsp>\d+)",
            )
            
            mstau = d_parse_result["mstau"]
            ctau = d_parse_result["ctau"]
            
            proclabel = rf"\\vtop{{\\hbox{{$m_{{\\tilde{{\\tau}}}}({mstau})$}}\\hbox{{$c\\tau_{{0}}({ctau})$}}}}"
        
        d_cutflow_output[prockey] = {
            "label" : str(proclabel),
            "cuts" : [],
        }
        scaleby_proc = eval(d_procinfo.get("scaleby", "1"))
        
        for icut, d_cutinfo in enumerate(d_config["cuts"]) :
            
            cut = d_cutinfo["key"]
            print(f"  [cut {cut}]")
            
            scaleby_cut = eval(d_cutinfo.get("scaleby", "1"))
            nevents_proc_cut = 0
            nevents_err_proc_cut = 0
            
            for sample in d_procinfo["samples"] :
                
                print(f"    [sample {sample}]")
                
                for d_erainfo in d_config["eras"] :
                    
                    era = d_erainfo["label"]
                    print(f"      [era {era}]")
                    
                    lumi = d_erainfo["lumi"]
                    scaleby_era = eval(d_erainfo.get("scaleby", "1"))
                    nevents_tot_era = functools.reduce(operator.getitem, [sample]+d_config["neventkey"].split("."), d_cutflow_input[era])
                    #print(d_erainfo)
                    
                    nevents_raw = functools.reduce(operator.getitem, [sample]+cut.split("."), d_cutflow_input[era])
                    efficiency = (nevents_raw * scaleby_proc * scaleby_era * scaleby_cut) / nevents_tot_era
                    nevents_proc_cut += efficiency * d_xsec[sample] * lumi
                    #nevents_proc_cut = nevents_proc_cut.item()
                    
                    # Binomial error
                    # Add in quadrature
                    nevents_err_proc_cut += (efficiency * (1 - efficiency) / nevents_tot_era) * (d_xsec[sample] * lumi)**2
                    #nevents_err_proc_cut += efficiency**2 * (1.0/nevents_raw + 1.0/nevents_tot_era)
                    
                    print(f"        [nevents_raw {nevents_raw}]")
            
            nevents_err_proc_cut = nevents_err_proc_cut**0.5
            
            print(f"  [proc {prockey}] [cut {cut}] [nevents {nevents_proc_cut:0.4g} +/- {nevents_err_proc_cut:0.4g}]")
            
            nevents_tot = nevents_proc_cut if icut == 0 else nevents_tot
            
            d_cutflow_output[prockey]["cuts"].append({
                "cut_label" : d_cutinfo["label"],
                "yield" : float(f"{nevents_proc_cut: 0.4g}"),
                "yield_unc_stat" : float(f"{nevents_err_proc_cut: 0.4g}"),
                "efficiency" : float(f"{nevents_proc_cut / nevents_tot: 0.4g}"),
                "efficiency_unc_stat" : float(f"{nevents_err_proc_cut/nevents_tot: 0.4g}"),
            })
    
    cmut.logger.info(f"Writing cutflows to: {outfname}")
    with open(outfname, "w", encoding = "utf-8") as fopen:
        
        yaml.dump(d_cutflow_output, fopen)


if __name__ == "__main__" :
    main()