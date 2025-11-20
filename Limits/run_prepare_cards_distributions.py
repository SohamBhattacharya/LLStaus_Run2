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
            #"2016_preVFP.2016_postVFP",
            "2017",
            "2018"
        ]
    )
    
    parser.add_argument(
        "--outsuffix",
        help = "Will append \"_<outsuffix>\" to the output directory",
        type = str,
        required = False,
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    d_distributions = {}
    
    # nbins = -1 will set nbins = len(rebin) - binstart
    
    d_distributions["jet2_pt"] = {
        "rebin": [
          0,
          30,
          50,
          100,
          300,
          500,
        ],
        "binstart": 1,
        "nbins": -1,
    }
    
    d_distributions["METpt"] = {
        "rebin": [
          0,
          120,
          150,
          200,
          250,
          350,
          500,
        ],
        "binstart": 1,
        "nbins": -1,
    }
    
    d_distributions["mt2_j1_j2_MET"] = {
        "rebin": [
          0,
          54,
          102,
          150,
          204,
          300,
        ],
        "binstart": 1,
        "nbins": -1,
    }
    
    d_config_files = {}
    
    d_config_files["2016_preVFP"] = [
        "configs/distributions/2016_preVFP/config_distribution_var_BRT2.yaml"
    ]
    
    d_config_files["2016_postVFP"] = [
        "configs/distributions/2016_postVFP/config_distribution_var_BRT2.yaml"
    ]
    
    d_config_files["2017"] = [
        "configs/distributions/2017/config_distribution_var_BRT2.yaml"
    ]
    
    d_config_files["2018"] = [
        "configs/distributions/2018/config_distribution_var_BRT2.yaml"
    ]
    
    l_config_files = []
    
    for era in args.eras :
        
        l_config_files.extend(d_config_files[era])
    
    username = getpass.getuser()
    tmp_dir = f"/tmp/{username}"
    
    outdir = "results/distributions"
    
    outdir = f"{outdir}_{args.outsuffix}" if len(args.outsuffix) else outdir
    os.system(f"mkdir -p {outdir}")
    
    for varname, dist_info in d_distributions.items() :
        
        cmut.logger.info(f"Starting to process [{varname}]")
        
        l_config_files_tmp = []
        
        for cfg_file in l_config_files:
            
            file_content = None
            
            with open(cfg_file) as fin :
                
                file_content = fin.read()
            
            formatinfo = {
                "varname": varname,
                "rebin": dist_info["rebin"],
                "binstart": dist_info["binstart"],
                "nbins": dist_info["nbins"] if dist_info["nbins"] > 0 else len(dist_info["rebin"]) - dist_info["binstart"],
            }
            
            file_content = file_content.format(**formatinfo)
            
            fname, fext = os.path.splitext(os.path.basename(cfg_file))
            os.system(f"mkdir -p {tmp_dir}")
            #fout_name = f"{tmp_dir}/{fname}_{varname}_{time.time_ns()}{fext}"
            fout_name = f"{tmp_dir}/{fname}_{varname}{fext}"
            fout_name = f"{tmp_dir}/{os.path.splitext(cfg_file)[0].replace('/', '_')}_{varname}{fext}"
            cmut.logger.info(f"Writing {cfg_file} --> {fout_name}")
            
            with open(fout_name, "w") as fout :
                
                fout.write(file_content)
                l_config_files_tmp.append(fout_name)
        
        #if (args.era != "2016_preVFP.2016_postVFP") :
        #    print(args.era)
        #    exit(1)
        
        cmd = " ".join([
            "set -x; python3 prepare_cards.py",
            f"--configs {' '.join(l_config_files_tmp)}",
            #f"--outdir {outdir}",
            "--yields_uncs",
            "--yields_uncs_sigs configs/distributions/sig_list_for-distributions.txt",
            "--nocards",
            "--combpars era",
            f"--outdir {outdir}",
        ])
        
        #cmd += f" --eracombos {args.era.replace('.', ',')}" * int(args.era == "2016_preVFP.2016_postVFP")
        
        cmd_retval = os.system(cmd)
        
        if (not cmd_retval) :
            
            cmut.logger.info(f"Processed [{varname}]")
        
        else :
            
            cmut.logger.info(f"Failure processing [{varname}]. Exiting...")
            exit(cmd_retval)
        
        #print("="*100)
        #print("\n\n")

if (__name__ == "__main__") :
    
    main()
