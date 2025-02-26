#!/usr/bin/env python3

import argparse
import os


def exec_cmds(cmds) :
    
    for cmd in cmds :
        
        retval = os.system(f"set -x; {cmd}")
        
        if (retval) :
            
            exit(retval)


def main() :
    
    # Argument parser
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        "--indir",
        help = "Input directory; can have wild cards, e.g.\"dirA/dirB*\"",
        type = str,
        required = True,
    )
    
    #parser.add_argument(
    #    "--outdir",
    #    help = "Output directory",
    #    type = str,
    #    required = True,
    #)
    
    parser.add_argument(
        "--wspace",
        help = "Create workspace",
        action = "store_true",
    )
    
    parser.add_argument(
        "--fit",
        help = "Run fit",
        action = "store_true",
    )
    
    parser.add_argument(
        "--collfits",
        help = "Collect fit results",
        action = "store_true",
    )
    
    parser.add_argument(
        "--scan",
        help = "Run NLL scan",
        action = "store_true",
    )
    
    parser.add_argument(
        "--postfit",
        help = "Run postfit",
        action = "store_true",
    )
    
    parser.add_argument(
        "--impacts",
        help = "Run impacts",
        action = "store_true",
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    indir = os.path.dirname(args.indir)
    
    cmds = []
    
    if (args.wspace) :
        
        cmd = (
            "nice -n 10 combineTool.py"
            " -M T2W"
            f" -i {args.indir}/card_ana_DisTauSF_mss_ZMT_*.txt"
            " -o workspace.root"
            " -m 90"
            " -P TauFW.Fitter.models.TagAndProbeModel:tagAndProbe"
            " --PO verbose=2"
            " --PO pass=pass"
            " --PO fail=fail"
            " --parallel 15"
        )
        
        cmds.append(cmd)
    
    if (args.fit) :
        
        cmd = (
            "nice -n 10 combineTool.py"
            " -M FitDiagnostics"
            f" -d {args.indir}/workspace.root"
            " --redefineSignalPOIs SF"
            " --setParameterRanges SF=0,3"
            " --saveShapes"
            " --saveWithUncertainties"
            " --saveOverallShapes"
            " --saveNormalizations"
            " --numToysForShapes 200"
            " --forceRecreateNLL"
            " --skipBOnlyFit"
            " --cminDefaultMinimizerStrategy 0"
            " --robustFit 1"
            " --plots"
            " --there"
            " -n \".DisTauSF\""
            " --parallel 15"
        )
        
        cmds.append(cmd)
    
    if (args.collfits) :
        
        cmd = (
            "nice -n 10 combineTool.py"
            " -M CollectLimits"
            f" {args.indir}/higgsCombine.DisTauSF.FitDiagnostics.mH120.root"
            " --use-dirs"
            f" -o {indir}/fit_result.json"
        )
        
        cmds.append(cmd)
    
    if (args.postfit) :
        
        cmd = (
            "./run_postfit_DisTauSF.sh"
            f" \"{args.indir}\""
        )
        
        cmds.append(cmd)
    
    if (args.scan) :
        
        wkdir_scan = "scan"
        cmd = f"for dir in $(find {args.indir} -mindepth 0 -maxdepth 0 -type d); do mkdir -pv $dir/{wkdir_scan}; cp -v $dir/workspace.root $dir/{wkdir_scan}/; done"
        cmds.append(cmd)
        
        cmd = (
            "nice -n 10 combineTool.py"
            " -M MultiDimFit"
            f" -d {args.indir}/{wkdir_scan}/workspace.root"
            " --saveSpecifiedFunc SF_fail"
            " --setParameterRanges SF=0,3"
            " --robustFit 1"
            " --forceRecreateNLL"
            " --algo grid"
            " --points 300"
            " --cminDefaultMinimizerStrategy 0"
            " --there"
            " -n \".DisTauSF\""
            " --parallel 15"
        )
        
        cmds.append(cmd)
    
    if (args.impacts) :
        
        cmd = (
            "./run_impacts_DisTauSF.sh"
            f" \"{args.indir}\""
        )
        
        cmds.append(cmd)
    
    
    exec_cmds(cmds)
    
    return 0


if (__name__ == "__main__") :
    
    main()
