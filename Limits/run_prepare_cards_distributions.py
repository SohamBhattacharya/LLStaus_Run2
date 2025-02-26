#!/usr/bin/env python3

import argparse
import getpass
import os
import time
import utils.commonutils as cmut

def main() :
    
    # Argument parser
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        "--eras",
        help = "Eras",
        type = str,
        nargs = "+",
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
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    d_distributions = {}
    
    d_distributions["METpt"] = {
        "rebin": [
          0,
          120,
          150,
          250,
          500,
        ],
        "binstart": 2,
        "nbins": 3,
    }
    
    d_config_files = {}
    
    d_config_files["2016_preVFP"] = [
        "configs/distributions/2016_preVFP/config_distrubution_var.yaml"
    ]
    
    d_config_files["2016_postVFP"] = [
        "configs/distributions/2016_postVFP/config_distrubution_var.yaml"
    ]
    
    d_config_files["2017"] = [
        "configs/distributions/2017/config_distrubution_var.yaml"
    ]
    
    d_config_files["2018"] = [
        "configs/distributions/2018/config_datacard_var.yaml"
    ]
    
    l_config_files = []
    
    for era in args.eras :
        
        l_config_files.append(d_config_files[era])
    
    username = getpass.getuser()
    tmp_dir = f"/tmp/{username}"
    
    outdir = f"results/distributions"
    
    outdir = f"{outdir}_{args.outsuffix}" if len(args.outsuffix) else outdir
    os.system(f"mkdir -p {outdir}")
    
    for dist_name, dist_info in d_distributions.items() :
        
        nametag = f"{key_wp}_{key_dxy}"
        
        cmut.logger.info(f"Preparing cards for [wp {key_wp}], [dxy {key_dxy}], [rebin {rebin}], [binstart {binstart}], [nbins {nbins}]")
        
        config_files_wp = []
        
        for cfg_file in l_config_files:
            
            file_content = None
            
            with open(cfg_file) as fin :
                
                file_content = fin.read()
            
            formatinfo = {
                "dist_name": dist_name,
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
