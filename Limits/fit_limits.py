#!/usr/bin/env python3

import argparse
import glob
import os

import utils.commonutils as cmut


def exec_cmds(cmds) :
    
    for cmd in cmds :
        
        print(cmd)
        retval = os.system(cmd)
        
        if (retval) :
            
            exit(retval)


def main() :
    
    # Argument parser
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        "--indir",
        help = "Input directory; can have wildcards, e.g.\"dirA/dirB*\"",
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
        "--limits",
        help = "Run limit calculation",
        action = "store_true",
    )
    
    parser.add_argument(
        "--fitdiag",
        help = "Run fit diagnostics",
        action = "store_true",
    )
    
    parser.add_argument(
        "--colllimits",
        help = "Collect limit results",
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
    
    parser.add_argument(
        "--plotcorr",
        help = "Plot bin-to-bin and nuisance correlations",
        action = "store_true",
    )
    
    parser.add_argument(
        "--signal",
        help = "Will pass --expectSignal <signal> to combine",
        type = str,
        required = False,
        choices = ["0", "1"],
        default = "",
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    nparallel = 20
    
    signal_suffix = f"_expectSignal{args.signal}" if len(args.signal) else ""
    signal_args = f"-t -1 --expectSignal {args.signal}" if len(args.signal) else ""
    
    indir = os.path.dirname(args.indir)
    
    cmds = []
    
    # Note: ordering is important!
    # For example, cannot workspace production must come first
    
    if (args.wspace) :
        
        cmd = (
            "nice -n 10 combineTool.py"
            " -M T2W"
            f" -i {args.indir}/card_*.txt"
            " -m 90"
            " -o workspace.root"
            " --PO verbose=2"
            f" --parallel {nparallel}"
        )
        
        cmds.append(cmd)
    
    if (args.limits) :
        
        wkdir_limits = "limits"
        cmd = f"for dir in $(find {args.indir} -mindepth 0 -maxdepth 0 -type d); do mkdir -pv $dir/{wkdir_limits}; cp -v $dir/workspace.root $dir/{wkdir_limits}/; done"
        cmds.append(cmd)
        
        cmd = (
            "nice -n 10 combineTool.py"
            " -M AsymptoticLimits"
            f" -d {args.indir}/{wkdir_limits}/workspace.root"
            " --there"
            f" --parallel {nparallel}"
        )
        
        cmds.append(cmd)
    
    if (args.colllimits) :
        
        #cmd = (
        #    "nice -n 10 combineTool.py"
        #    " -M CollectLimits"
        #    f" {args.indir}/limits/higgsCombine.Test.AsymptoticLimits.mH120.root"
        #    #" --use-dirs"
        #    f" -o {indir}/limits/limits.json"
        #)
        
        wkdir_limits = "limits"
        
        cmd = (
            f"for dir in $(find {args.indir} -mindepth 0 -maxdepth 0 -type d); do \n"
            f"  nice -n 10 combineTool.py -M CollectLimits $dir/{wkdir_limits}/higgsCombine.Test.AsymptoticLimits.mH120.root -o $dir/{wkdir_limits}/limits.json & \n"
            f"  while [[ $(jobs -rp | wc -l) -ge {nparallel} ]]; do \n"
            "       sleep 1 \n"
            "   done \n"
            " done \n"
            " wait $(jobs -rp) \n"
            " echo"
        )
        
        cmds.append(cmd)
    
    if (args.postfit) :
        
        cmd = (
            "./run_postfit_limits.sh"
            f" \"{args.indir}\""
        )
        
        cmds.append(cmd)
    
    if (args.fitdiag) :
        
        wkdir_fitdiag = f"fit-diagnostics{signal_suffix}"
        cmd = f"for dir in $(find {args.indir} -mindepth 0 -maxdepth 0 -type d); do mkdir -pv $dir/{wkdir_fitdiag}; cp -v $dir/workspace.root $dir/{wkdir_fitdiag}/; done"
        cmds.append(cmd)
        
        cmd = (
            "nice -n 10 combineTool.py"
            " -M FitDiagnostics"
            f" -d {args.indir}/{wkdir_fitdiag}/workspace.root"
            " --saveShapes"
            " --saveWithUncertainties"
            " --saveOverallShapes"
            " --saveNormalizations"
            " --numToysForShapes 200"
            " --forceRecreateNLL"
            " --robustFit 1"
            " --setRobustFitAlgo Minuit2"
            " --maxFailedSteps 20"
            " --rMin -100"
            f" {signal_args}"
            " --plots"
            " --there"
            #" -n \".limits\""
            f" --parallel {nparallel}"
        )
        
        cmds.append(cmd)
    
    if (args.scan) :
        
        wkdir_scan = f"scan{signal_suffix}"
        cmd = f"for dir in $(find {args.indir} -mindepth 0 -maxdepth 0 -type d); do mkdir -pv $dir/{wkdir_scan}; cp -v $dir/workspace.root $dir/{wkdir_scan}/; done"
        cmds.append(cmd)
        
        cmd = (
            "nice -n 10 combineTool.py"
            " -M MultiDimFit"
            f" -d {args.indir}/{wkdir_scan}/workspace.root"
            #" --setParameterRanges r=-3,3"
            " --robustFit 1"
            " --setRobustFitAlgo Minuit2"
            " --maxFailedSteps 20"
            " --setRobustFitStrategy=2"
            " --setRobustFitTolerance=0.01"
            " --rMin -3"
            " --rMax 3"
            " --forceRecreateNLL"
            #" --cminDefaultMinimizerStrategy 0"
            " --algo grid"
            " --points 500"
            f" {signal_args}"
            " --there"
            #" -n \".scan\""
            f" --parallel {nparallel}"
        )
        
        cmds.append(cmd)
    
    if (args.impacts) :
        
        cmd = (
            "./run_impacts_limits.sh"
            f" \"{args.indir}\""
            f" {args.signal}"
        )
        
        cmds.append(cmd)
    
    
    exec_cmds(cmds)
    
    
    if (args.plotcorr) :
        
        wkdir = f"{args.indir}/{wkdir_fitdiag}" if "wkdir_fitdiag" in locals() else args.indir
        
        l_fnames = glob.glob(f"{wkdir}/fitDiagnostics.Test.root")
        
        l_histnames = [
            "shapes_prefit/overall_total_covar",
            "shapes_fit_b/overall_total_covar",
            "shapes_fit_s/overall_total_covar",
        ]
        
        d_info = {
            "shapes_prefit/overall_total_covar": {
                "outfilename": "bin-to-bin_correlation_prefit",
                "title": "Bin-to-bin correlation (prefit)",
                "cov_to_corr": True,
            },
            
            "shapes_fit_b/overall_total_covar": {
                "outfilename": "bin-to-bin_correlation_fit_b",
                "title": "Bin-to-bin correlation (background only fit)",
                "cov_to_corr": True,
            },
            
            "shapes_fit_s/overall_total_covar": {
                "outfilename": "bin-to-bin_correlation_fit_sb",
                "title": "Bin-to-bin correlation (signal+background fit)",
                "cov_to_corr": True,
            },
            
            "covariance_fit_b": {
                "outfilename": "nuisance_correlation_fit_b",
                "title": "Nuisance correlation (background only fit)",
                "cov_to_corr": False,
            },
            
            "covariance_fit_s": {
                "outfilename": "nuisance_correlation_fit_sb",
                "title": "Nuisance correlation (signal+background fit)",
                "cov_to_corr": False,
            }
        }
        
        for fname in l_fnames :
            
            for histname, info in d_info.items() :
                
                cmut.logger.info(f"Plotting correlation from [hist {histname}] [file {fname}]")
                
                cmut.plot_fitdiagnostics_correlation(
                    rootfilename = fname,
                    outfilename = f"{os.path.dirname(fname)}/{info['outfilename']}.pdf",
                    histname = histname,
                    l_binlabel_fmt = [
                        "{label}.replace('BRT2_', '')",
                        "{label}.replace('_0', '')",
                    ],
                    title = info["title"],
                    cov_to_corr = info["cov_to_corr"]
                )
    
    return 0


if (__name__ == "__main__") :
    
    main()

