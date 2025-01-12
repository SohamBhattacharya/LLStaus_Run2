#!/usr/bin/env python3

import argparse
import getpass
import os
import time
import utils.commonutils as cmut

def main() :
    
    wps = [
        "p9900",
        "p9800",
        "p9700",
        "p9500",
        "p9000",
        "p8000",
    ]
    
    # ["name", "binning", binstart]
    binnings = {
        "dxy-gt0p00": {
            "rebin": "[0, 2]",
            "binstart": 1,
            "nbins": 1,
        },
        "dxy-gt0p01": {
            "rebin": "[0, 0.01, 2]",
            "binstart": 2,
            "nbins": 1,
        },
        "dxy-gt0p02": {
            "rebin": "[0, 0.02, 2]",
            "binstart": 2,
            "nbins": 1,
        },
        "dxy-gt0p03": {
            "rebin": "[0, 0.03, 2]",
            "binstart": 2,
            "nbins": 1,
        },
        "dxy-gt0p04": {
            "rebin": "[0, 0.04, 2]",
            "binstart": 2,
            "nbins": 1,
        },
    }
    
    config_files = [
        "configs/config_datacard_mutau_pass_dxy.yaml",
        "configs/config_datacard_mutau_fail_dxy.yaml",
        "configs/config_datacard_mumu.yaml",
    ]
    
    username = getpass.getuser()
    tmp_dir = f"/tmp/{username}"
    
    for iwp, wp in enumerate(wps) :
        
        for binname, bininfo in binnings.items() :
            
            rebin = bininfo["rebin"]
            binstart = bininfo["binstart"]
            nbins = bininfo["nbins"]
            
            nametag = f"{wp}_{binname}"
            
            cmut.logger.info(f"Preparing cards for [wp {wp}], [rebin {rebin}], [binstart {binstart}], [nbins {nbins}]")
            
            config_files_wp = []
            
            for cfg_file in config_files:
                
                file_content = None
                
                with open(cfg_file) as fin :
                    
                    file_content = fin.read()
                
                formatinfo = {
                    "wp": wp,
                    "nametag": nametag,
                }
                
                formatinfo.update(bininfo)
                
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
                " --outdir tmp/test_DisTauSF"
            )
            
            cmd_retval = os.system(cmd)
            
            if (not cmd_retval) :
                
                cmut.logger.info(f"Prepared cards for wp {wp}")
            
            else :
                
                cmut.logger.info(f"Failure preparing cards for wp {wp}. Exiting...")
                exit(cmd_retval)
            
            print("="*100)
            print("\n\n")

if (__name__ == "__main__") :
    
    main()
