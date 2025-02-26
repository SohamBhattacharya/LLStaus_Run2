#!/usr/bin/env python3

import argparse
import getpass
import os
import time

from traitlets import default
import utils.commonutils as cmut

def main() :
    
    # Argument parser
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        "--era",
        help = "Era",
        type = str,
        required = True,
        choices = [
            "2016_preVFP",
            "2016_postVFP",
            "2016_preVFP.2016_postVFP",
            "2017",
            "2018"
        ]
    )
    
    parser.add_argument(
        "--outsuffix",
        help = "Will append \"_outsuffix\" to the output directory",
        type = str,
        required = False,
        default = ""
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # label: dir name
    d_wp = {
        ##"wp-p00": "all",
        
        "wp-p70": "p7000",
        "wp-p80": "p8000",
        "wp-p90": "p9000",
        "wp-p95": "p9500",
        "wp-p97": "p9700",
        "wp-p98": "p9800",
        "wp-p99": "p9900",
    }
    
    d_dxy = {
        "dxy-gt-0p00"    : "all",
        "dxy-gt-0p01"    : "dxy01",
        "dxy-gt-0p02"    : "dxy02",
        "dxy-gt-0p03"    : "dxy03",
        "dxy-gt-0p05"    : "dxy05",
        "dxy-gt-0p07"    : "dxy07",
        "dxy-gt-0p09"    : "dxy09",
    }
    
    d_config_files = {}
    
    d_config_files["2016_preVFP"] = [
        "configs/DisTauSF/mutau_mass/2016_preVFP/config_datacard_mutau_pass_mass.yaml",
        "configs/DisTauSF/mutau_mass/2016_preVFP/config_datacard_mutau_fail_mass.yaml",
    ]
    
    d_config_files["2016_postVFP"] = [
        "configs/DisTauSF/mutau_mass/2016_postVFP/config_datacard_mutau_pass_mass.yaml",
        "configs/DisTauSF/mutau_mass/2016_postVFP/config_datacard_mutau_fail_mass.yaml",
    ]
    
    d_config_files["2016_preVFP.2016_postVFP"] = d_config_files["2016_preVFP"] + d_config_files["2016_postVFP"]
    
    d_config_files["2017"] = [
        "configs/DisTauSF/mutau_mass/2017/config_datacard_mutau_pass_mass.yaml",
        "configs/DisTauSF/mutau_mass/2017/config_datacard_mutau_fail_mass.yaml",
        
        #"configs/DisTauSF/mutau_mass/2017/config_datacard_mutau_pass_mass_w-distau-WPL.yaml",
        #"configs/DisTauSF/mutau_mass/2017/config_datacard_mutau_fail_mass_w-distau-WPL.yaml",
    ]
    
    d_config_files["2018"] = [
        #"configs/DisTauSF/mutau_mass/2018/config_datacard_mutau_pass_mass.yaml",
        #"configs/DisTauSF/mutau_mass/2018/config_datacard_mutau_fail_mass.yaml",
        
        #"configs/DisTauSF/mumu/2018/config_datacard_mumu.yaml",
        
        "configs/DisTauSF/mutau_mass/2018/DeepTau-Loose/config_datacard_mutau_pass_mass_DeepTau-Loose.yaml",
        "configs/DisTauSF/mutau_mass/2018/DeepTau-Loose/config_datacard_mutau_fail_mass_DeepTau-Loose.yaml",
    ]
    
    l_config_files = d_config_files[args.era]
    
    #outsuffix = ""
    #outsuffix = "mumu-scaled-1em2"
    #outsuffix = "mumu-scaled-1em3"
    #outsuffix = "mumu-scaled-1em4"
    #outsuffix = "no-dy-rateparam"
    
    username = getpass.getuser()
    tmp_dir = f"/tmp/{username}"
    
    # Binning
    
    rebin = [60, 80]
    binstart = 1
    nbins = 1
    
    #rebin = [60, 70, 80]
    #binstart = 1
    #nbins = 2
    
    #rebin = [60, 66, 72, 80]
    #binstart = 1
    #nbins = 3
    
    #outdir = f"tmp/test_DisTauSF_mass_nbins{nbins}-{'-'.join([str(_ele) for _ele in rebin])}"
    outdir = f"results/DisTauSF/DisTauSF_mass_nbins{nbins}-{'-'.join([str(_ele) for _ele in rebin])}"
    
    outdir = f"{outdir}_{args.outsuffix}" if args.outsuffix else outdir
    os.system(f"mkdir -p {outdir}")
    
    for key_wp, val_wp in d_wp.items() :
        
        for key_dxy, val_dxy in d_dxy.items() :
            
            nametag = f"{key_wp}_{key_dxy}"
            
            cmut.logger.info(f"Preparing cards for [wp {key_wp}], [dxy {key_dxy}], [rebin {rebin}], [binstart {binstart}], [nbins {nbins}]")
            
            config_files_wp = []
            
            for cfg_file in l_config_files:
                
                file_content = None
                
                with open(cfg_file) as fin :
                    
                    file_content = fin.read()
                
                formatinfo = {
                    "wp": val_wp,
                    "dxy": val_dxy,
                    "nametag": nametag,
                    "rebin": rebin,
                    "binstart": binstart,
                    "nbins": nbins,
                }
                
                file_content = file_content.format(**formatinfo)
                
                fname, fext = os.path.splitext(os.path.basename(cfg_file))
                os.system(f"mkdir -p {tmp_dir}")
                fout_name = f"{tmp_dir}/{fname}_{key_wp}_{key_dxy}_{time.time_ns()}{fext}"
                cmut.logger.info(f"Writing {cfg_file} --> {fout_name}")
                
                with open(fout_name, "w") as fout :
                    
                    fout.write(file_content)
                    config_files_wp.append(fout_name)
            
            #if (args.era != "2016_preVFP.2016_postVFP") :
            #    print(args.era)
            #    exit(1)
            
            cmd = (
                "set -x; python3 prepare_cards.py"
                f" --configs {' '.join(config_files_wp)}"
                f" --outdir {outdir}"
                f" --combpars channel"
                #" --yields_uncs"
                #f" --chcombos mutau_pass,mutau_fail"
            )
            
            cmd += f" --eracombos {args.era.replace('.', ',')}" * int(args.era == "2016_preVFP.2016_postVFP")
            
            cmd_retval = os.system(cmd)
            
            if (not cmd_retval) :
                
                cmut.logger.info(f"Prepared cards for [wp {key_wp}], [dxy {key_dxy}] in {outdir}")
            
            else :
                
                cmut.logger.info(f"Failure preparing cards for [wp {key_wp}], [dxy {key_dxy}]. Exiting...")
                exit(cmd_retval)
            
            #print("="*100)
            #print("\n\n")

if (__name__ == "__main__") :
    
    main()
