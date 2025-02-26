#!/usr/bin/env python3

import argparse
import itertools
import os


ENABLE_OPT = "enable"
DISABLE_OPT = "disable"


def main() :
    
    # Argument parser
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        "--indir",
        help = "Input directory; can have wildcards, e.g.\"dirA/dirB*\"",
        type = str,
        required = True,
    )
    
    parser.add_argument(
        "--mode",
        help = "Enable or disable each bin",
        type = str,
        required = True,
        choices = [ENABLE_OPT, DISABLE_OPT]
    )
    
    parser.add_argument(
        "--wspace",
        help = "Create workspace",
        action = "store_true",
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    
    l_eras = [
        "2016_preVFP",
        "2016_postVFP",
        "2017",
        "2018",
    ]

    binnames = [
        "bin3_BRT2",
        "bin4_BRT2",
        "bin5_BRT2",
        "bin6_BRT2",
        "bin7_BRT2",
        "bin8_BRT2",
        "bin9_BRT2",
        "bin10_BRT2",
    ]
    
    if args.wspace :
        
        cmd_retval = os.system(f"./fit_limits.py --wspace --indir \"{args.indir}\"")
        
        if(cmd_retval) :
            
            exit(cmd_retval)
    
    for binname in binnames :
        
        suffix = None
        
        if args.mode == ENABLE_OPT :
            
            binnames_mask = [[f"mask_{_bin}_{_era}=1" for _era in l_eras] for _bin in binnames if binname != _bin]
            suffix = f"enable_{binname}"
        
        elif args.mode == DISABLE_OPT :
            
            binnames_mask = [[f"mask_{_bin}_{_era}=1" for _era in l_eras] for _bin in binnames if binname == _bin]
            suffix = f"disable_{binname}"
        
        binnames_mask = list(itertools.chain.from_iterable(binnames_mask))
        mask_str = ",".join(binnames_mask)
        
        assert (suffix is not None)
        
        cmd = " ".join([
            "./fit_limits.py",
            #"--help",
            #"--wspace",
            "--limits",
            "--colllimits",
            f"--combineargs \"--setParameters {mask_str}\"",
            f"--suffix {suffix}",
            f"--indir \"{args.indir}\"",
        ])
        
        print(cmd)
        
        cmd_retval = os.system(cmd)
        
        if(cmd_retval) :
            
            exit(cmd_retval)
    
    return 0


if __name__ == "__main__" :
    
    main()