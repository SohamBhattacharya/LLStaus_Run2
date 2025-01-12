#!/usr/bin/env python3

import argparse
import os


def exec_cmds(cmds) :
    
    for cmd in cmds :
        
        retval = os.system(cmd)
        
        if (retval) :
            
            exit(retval)


def main() :
    
    # Argument parser
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        "--indir",
        help = "Input directory; can have wild cards, e.g.""dirA/dirB*\"",
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
            "combineTool.py"
            " -M T2W"
            f" -i {args.indir}/card_DisTauSF_ZMT_*.txt"
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
            "combineTool.py"
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
            " --plots"
            " --there"
            " -n \".DisTauSF\""
            " --parallel 15"
        )
        
        cmds.append(cmd)
    
    if (args.collfits) :
        
        cmd = (
            "combineTool.py"
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
        
        cmd = (
            "combineTool.py"
            " -M MultiDimFit"
            f" -d {args.indir}/workspace.root"
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

#!/bin/bash

# Exit if a command fails
#set -Eeu -o pipefail

# Prints command before executing
#set -o xtrace


#DIR="$1"
#DIR=tmp/test_DisTauSF_mass_nbins1-60-80/DisTauSF/channels_all/eras_all/ZMT_wp-*_dxy-gt-*
#DIR=tmp/test_DisTauSF_mass_nbins2-60-70-80/DisTauSF/channels_all/eras_all/ZMT_wp-*_dxy-gt-*
#DIR=tmp/test_DisTauSF_mass_nbins3-60-66-72-80/DisTauSF/channels_all/eras_all/ZMT_wp-*_dxy-gt-*
#OUTDIR=tmp/test_DisTauSF_mass_nbins1-60-80
#OUTDIR=$(dirname "$DIR")

#echo $DIR
#echo $OUTDIR
#exit

#combineTool.py"
#-M T2W"
#-i $DIR/card_DisTauSF_ZMT_*.txt"
#-o workspace.root"
#-m 90"
#-P TauFW.Fitter.models.TagAndProbeModel:tagAndProbe"
#--PO verbose=2"
#--PO pass=pass"
#--PO fail=fail"
#--parallel 15


#combineTool.py"
#-M FitDiagnostics"
#-d $DIR/workspace.root"
#--redefineSignalPOIs SF"
#--setParameterRanges SF=0,3"
#--saveShapes"
#--saveWithUncertainties"
#--saveOverallShapes"
#--numToysForShapes 200"
#--forceRecreateNLL"
#--skipBOnlyFit"
#--cminDefaultMinimizerStrategy 0"
#--plots"
#--there"
#-n ".DisTauSF""
#--parallel 15


#combineTool.py"
#-M MultiDimFit"
#-d $DIR/workspace.root"
#--saveSpecifiedFunc SF_fail"
#--setParameterRanges SF=0,3"
#--robustFit 1"
#--forceRecreateNLL"
#--algo grid"
#--points 300"
#--cminDefaultMinimizerStrategy 0"
#--there"
#-n ".DisTauSF""
#--parallel 15

#--saveSpecifiedFunc SF_fail,rp_dy_norm_2018"
#--cl 0.6827"

#
#combineTool.py"
#-M CollectLimits"
#$DIR/higgsCombine.DisTauSF.FitDiagnostics.mH120.root"
#--use-dirs"
#-o $OUTDIR/fit_result.json