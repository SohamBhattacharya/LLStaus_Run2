#!/usr/bin/env python3

import argparse
import getpass
import os
import time
import utils.commonutils as cmut

def main() :
    
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
        "dxy-gt-0p01"    : "dxy0p1",
        "dxy-gt-0p02"    : "dxy0p2",
        "dxy-gt-0p03"    : "dxy0p3",
        "dxy-gt-0p05"    : "dxy0p5",
        "dxy-gt-0p07"    : "dxy0p7",
        "dxy-gt-0p09"    : "dxy0p9",
    }
    
    config_files = [
        "configs/config_datacard_mutau_pass_mass.yaml",
        "configs/config_datacard_mutau_fail_mass.yaml",
        "configs/config_datacard_mumu.yaml",
    ]
    
    #out_suffix = ""
    out_suffix = "mumu-scaled-1em2"
    #out_suffix = "mumu-scaled-1em3"
    #out_suffix = "mumu-scaled-1em4"
    #out_suffix = "no-dy-rateparam"
    
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
    
    outdir = f"tmp/test_DisTauSF_mass_nbins{nbins}-{'-'.join([str(_ele) for _ele in rebin])}"
    
    outdir = f"{outdir}_{out_suffix}" if len(out_suffix) else outdir
    
    for key_wp, val_wp in d_wp.items() :
        
        for key_dxy, val_dxy in d_dxy.items() :
            
            nametag = f"{key_wp}_{key_dxy}"
            
            cmut.logger.info(f"Preparing cards for [wp {key_wp}], [dxy {key_dxy}], [rebin {rebin}], [binstart {binstart}], [nbins {nbins}]")
            
            config_files_wp = []
            
            for cfg_file in config_files:
                
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
                fout_name = f"{tmp_dir}/{fname}_{time.time_ns()}{fext}"
                cmut.logger.info(f"Writing {cfg_file} --> {fout_name}")
                
                with open(fout_name, "w") as fout :
                    
                    fout.write(file_content)
                    config_files_wp.append(fout_name)
            
            cmd = (
                "python3 prepare_cards.py "
                f"--configs {' '.join(config_files_wp)}"
                f" --outdir {outdir}"
                f" --chcombos mutau_pass,mutau_fail"
            )
            
            cmd_retval = os.system(cmd)
            
            if (not cmd_retval) :
                
                cmut.logger.info(f"Prepared cards for [wp {key_wp}], [dxy {key_dxy}] in {outdir}")
            
            else :
                
                cmut.logger.info(f"Failure preparing cards for [wp {key_wp}], [dxy {key_dxy}]. Exiting...")
                exit(cmd_retval)
            
            print("="*100)
            print("\n\n")

if (__name__ == "__main__") :
    
    main()
