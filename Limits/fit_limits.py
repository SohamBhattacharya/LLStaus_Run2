#!/usr/bin/env python3

import argparse
import datetime
import glob
import itertools
import json
import numpy
import os
import time
import ROOT

import utils.commonutils as cmut


def main() :
    
    #timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    
    # Argument parser
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        "--indir",
        help = "Input directory; can have wildcards, e.g.\"dirA/dirB*\"",
        type = str,
        required = True,
    )
    
    parser.add_argument(
        "--suffix",
        help = "Will append \"_suffix\" to the combine output directories (only implemented for limits, fit diagnostics, and NLL scan)",
        type = str,
        required = False,
    )
    
    parser.add_argument(
        "--combineargs",
        help = (
            "Will pass this string as an argument to combine (only implemented for limits, fit diagnostics, and NLL scan); "
            "can contain multiple combine arguments in quotes: \"--arg1 --arg2 xyz\""
        ),
        type = str,
        required = False,
        default = ""
    )
    
    parser.add_argument(
        "--nparallel",
        help = "Number of parallel processes to run",
        type = int,
        required = False,
        default = 20
    )
    
    parser.add_argument(
        "--wspace",
        help = "Create workspace",
        action = "store_true",
    )
    
    parser.add_argument(
        "--limits",
        help = "Run asymptotic limit calculation",
        action = "store_true",
    )
    
    parser.add_argument(
        "--toylimits",
        help = "Run limit calculation with toys",
        action = "store_true",
    )
    
    #parser.add_argument(
    #    "--splitpoints",
    #    help = (
    #        "Run limit calculation with toys by generating a grid of test statistic distributions at various values of the signal strength. "
    #        "This will use the result of asymptotic limits to guess the range of signal strengths."
    #    ),
    #    action = "store_true",
    #)
    
    parser.add_argument(
        "--colllimits",
        help = "Collect limit results",
        action = "store_true",
    )
    
    parser.add_argument(
        "--haddtoylimits",
        help = "hadd the toy limit results from condor jobs",
        action = "store_true",
    )
    
    parser.add_argument(
        "--colltoylimits",
        help = (
            "Collect toy limit results. "
            "The --haddtoylimits argument must be passed if the files have not been hadded already."
            "The --timestamp argument must be passed to specify the condor job timestamp to read."
        ),
        action = "store_true",
    )
    
    parser.add_argument(
        "--fitdiag",
        help = "Run fit diagnostics",
        action = "store_true",
    )
    
    parser.add_argument(
        "--noplots",
        help = "Will not pass --plots to FitDiagnostics",
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
        "--gof",
        help = "Run goodness of fit tests",
        action = "store_true",
    )
    
    #parser.add_argument(
    #    "--collgof",
    #    help = "Collect goodness of fit test results",
    #    action = "store_true",
    #)
    
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
    
    parser.add_argument(
        "--timestamp",
        help = "Will use this timestamp (e.g. condor jobs for toys) instead of the default current timestamp. Usual format is YYYY-MM-DD_hh-mm-ss",
        type = str,
        required = False,
        default = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S"),
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    suffix = f"_{args.suffix}" if args.suffix else ""
    signal_suffix = f"_expectSignal{args.signal}" if len(args.signal) else ""
    signal_args = f"-t -1 --expectSignal {args.signal}" if len(args.signal) else ""
    
    wkdir_limits = f"limits{suffix}"
    wkdir_toylimits = f"toylimits{suffix}"
    wkdir_fitdiag = f"fit-diagnostics{suffix}{signal_suffix}"
    wkdir_scan = f"scan{suffix}{signal_suffix}"
    wkdir_gof = "gof"
    wkdir_plotcorr = f"{args.indir}/{wkdir_fitdiag}" if args.fitdiag else args.indir
    
    #indir = os.path.dirname(args.indir)
    
    cmds = []
    
    # Note: ordering is important!
    # For example, cannot workspace production must come first
    
    if (args.wspace) :
        
        cmd = (
            "nice -n 10 combineTool.py"
            " -M T2W"
            f" -i {args.indir}/card_*.txt"
            " --channel-masks"
            " -m 120"
            " -o workspace.root"
            " --PO verbose=2"
            f" --parallel {args.nparallel}"
            f" {args.combineargs}"
        )
        
        cmds.append(cmd)
    
    if (args.limits) :
        
        cmd = f"for dir in $(find {args.indir} -mindepth 0 -maxdepth 0 -type d); do mkdir -pv $dir/{wkdir_limits}; cp -v $dir/workspace.root $dir/{wkdir_limits}/; done"
        cmds.append(cmd)
        
        cmd = (
            "nice -n 10 combineTool.py"
            " -M AsymptoticLimits"
            f" -d {args.indir}/{wkdir_limits}/workspace.root"
            " --there"
            f" --parallel {args.nparallel}"
            f" {args.combineargs}"
        )
        
        cmds.append(cmd)
    
    
    l_toylimits_skipped = []
    
    if (args.toylimits) :
        
        d_limtypes = {
            "obs": "",
            "exp2p5": "--expectedFromGrid=0.025",
            "exp16": "--expectedFromGrid=0.16",
            "exp50": "--expectedFromGrid=0.5",
            "exp84": "--expectedFromGrid=0.84",
            "exp97p5": "--expectedFromGrid=0.975",
        }
        
        d_mappings_asymptotic = {
            "obs": 5,
            "exp0": 2,
            "exp-1": 1,
            "exp+1": 3,
            "exp-2": 0,
            "exp+2": 4,
        }
        
        l_indirs = cmut.natural_sort(glob.glob(f"{args.indir}*/"))
        d_limits_asymptotic = {}
        
        for indir in l_indirs :
            
            sig_info = cmut.parse_string_regex(
                s = indir,
                regexp = "SMS-TStauStau_MStau-(?P<mstau>\w+)_ctau-(?P<ctau>\w+)_mLSP-(?P<mlsp>\d+)",
            )
            
            sig_key = (sig_info["mstau"], sig_info["ctau"])
            
            try :
                rdf_asymptotic = ROOT.RDataFrame("limit", f"{indir}/limits/higgsCombine.Test.AsymptoticLimits.mH120.root")
                arr_limits_asymptotic = list(rdf_asymptotic.AsNumpy(["limit"])["limit"])
                
                d_limits_asymptotic[sig_key] = [_val for _val in arr_limits_asymptotic if numpy.isfinite(_val)]
                
                assert len(d_limits_asymptotic[sig_key]) > 0, f"No valid asymptotic limits found in: {arr_limits_asymptotic}"
                #assert len(d_limits_asymptotic[sig_key]), f"No asymptotic limits found for {sig_key}"
            
            except Exception as excp :
                
                l_toylimits_skipped.append((sig_key, indir, excp))
                cmut.logger.warning(f"Did not find asymptotic limits for {sig_key}; skipped. Total skipped: {len(l_toylimits_skipped)}")
                continue
            
            wkdir_tmp = f"{wkdir_toylimits}"
            condor_dir = f"{indir}/{wkdir_tmp}/condor/{args.timestamp}"
            
            cmd = " ".join([
                f"mkdir -pv {condor_dir};",
                f"cp -v {indir}/workspace.root {condor_dir};",
            ])
            cmds.append(cmd)
            
            #lim_min = min(arr_limits_asymptotic)
            lim_min = 0
            lim_max = round(2 * max(arr_limits_asymptotic), 8)
            lim_step = 0.01 * (lim_max - lim_min) # 1% of the range
            njobs = 10
            nmerge = max(int((lim_max - lim_min) / lim_step / njobs), 1)
            
            cmd = " ".join([
                f"pushd {condor_dir};"
                "combineTool.py",
                "-M HybridNew --LHCmode LHC-limits",
                "--there",
                "-d workspace.root",
                "--job-dir `pwd`",
                "--clsAcc 0",
                f"--singlePoint {lim_min}:{lim_max}:{lim_step}",
                f"--rMax {numpy.ceil(lim_max)}",
                f"--merge {nmerge}",
                "-T 20000",
                "--fullBToys",
                "-s -1",
                "--saveToys",
                "--saveHybridResult",
                "-m 120",
                "--task-name toylimits",
                "--job-mode condor",
                #"--sub-opts='+RequestRuntime=86400\nRequestMemory=6GB'",
                "--sub-opts='+RequestRuntime=259200'",
                args.combineargs,
                #"--dry-run",
                "; popd",
            ])
            
            cmds.append(cmd)
        
        #for limtype, limarg in d_limtypes.items() :
        #    
        #    wkdir_tmp = f"{wkdir_toylimits}_{limtype}"
        #    condor_dir = f"$dir/{wkdir_tmp}/condor/{args.timestamp}"
        #    
        #    cmd = " ".join([
        #        f"for dir in $(find {args.indir} -mindepth 0 -maxdepth 0 -type d | sort -V); do",
        #        #f"mkdir -pv $dir/{wkdir_tmp};",
        #        #f"cp -v $dir/workspace.root $dir/{wkdir_tmp}/;",
        #        f"mkdir -pv {condor_dir};",
        #        f"cp -v $dir/workspace.root {condor_dir};",
        #        "done",
        #    ])
        #    cmds.append(cmd)
        #    
        #    cmd = " ".join([
        #        f"for dir in $(find {args.indir} -mindepth 0 -maxdepth 0 -type d | sort -V); do",
        #        f"pushd {condor_dir};"
        #        "combineTool.py",
        #        "-M HybridNew --LHCmode LHC-limits",
        #        "--saveHybridResult",
        #        "--there",
        #        "-d workspace.root",
        #        "--job-dir `pwd`",
        #        limarg,
        #        #"--clsAcc 0.05",
        #        f"--task-name toylimits_{limtype}",
        #        args.combineargs,
        #        "--job-mode condor",
        #        #"--sub-opts='getenv=True\nrequest_memory=10GB\n+RequestRuntime=7200'",
        #        "--sub-opts='+RequestRuntime=86400\nRequestMemory=6GB'",
        #        #"--dry-run",
        #        "; popd",
        #        "; done",
        #    ])
        #    
        #    cmds.append(cmd)
    
    if (args.colllimits) :
        
        #cmd = (
        #    "nice -n 10 combineTool.py"
        #    " -M CollectLimits"
        #    f" {args.indir}/limits/higgsCombine.Test.AsymptoticLimits.mH120.root"
        #    #" --use-dirs"
        #    f" -o {indir}/limits/limits.json"
        #)
        
        cmd = (
            f"for dir in $(find {args.indir} -mindepth 0 -maxdepth 0 -type d | sort -V); do \n"
            f"  nice -n 10 combineTool.py -M CollectLimits $dir/{wkdir_limits}/higgsCombine.Test.AsymptoticLimits.mH120.root -o $dir/{wkdir_limits}/limits.json & \n"
            f"  while [[ $(jobs -rp | wc -l) -ge {args.nparallel} ]]; do \n"
            "       sleep 1 \n"
            "   done \n"
            " done \n"
            " wait $(jobs -rp) \n"
            " echo"
        )
        
        cmds.append(cmd)
    
    if (args.haddtoylimits) :
        
        cmd = (
            "./run_toylimits_from-grid.sh"
            f" \"{args.indir}/{wkdir_toylimits}/condor/{args.timestamp}\""
        )
        
        cmds.append(cmd)
    
    if (args.fitdiag) :
        
        cmd = f"for dir in $(find {args.indir} -mindepth 0 -maxdepth 0 -type d); do mkdir -pv $dir/{wkdir_fitdiag}; cp -v $dir/workspace.root $dir/{wkdir_fitdiag}/; done"
        cmds.append(cmd)
        
        cmd = " ".join([
            "nice -n 10 combineTool.py",
            "-M FitDiagnostics",
            f"-d {args.indir}/{wkdir_fitdiag}/workspace.root",
            "--saveShapes",
            "--saveWithUncertainties",
            "--saveOverallShapes",
            "--saveNormalizations",
            "--numToysForShapes 200",
            "--forceRecreateNLL",
            "--robustFit 1",
            "--setRobustFitAlgo Minuit2",
            "--maxFailedSteps 100",
            "--rMin -100",
            f" {signal_args}",
            "--plots"*(not args.noplots),
            "--there",
            #" -n \".limits\"",
            f"--parallel {args.nparallel}",
            f"{args.combineargs}",
        ])
        
        cmds.append(cmd)
    
    if (args.postfit) :
        
        cmd = (
            "./run_postfit_limits.sh"
            f" \"{args.indir}/{wkdir_fitdiag}\""
        )
        
        cmds.append(cmd)
    
    if (args.scan) :
        
        cmd = f"for dir in $(find {args.indir} -mindepth 0 -maxdepth 0 -type d); do mkdir -pv $dir/{wkdir_scan}; cp -v $dir/workspace.root $dir/{wkdir_scan}/; done"
        cmds.append(cmd)
        
        cmd = (
            "nice -n 10 combineTool.py"
            " -M MultiDimFit"
            f" -d {args.indir}/{wkdir_scan}/workspace.root"
            #" --setParameterRanges r=-3,3"
            " --robustFit 1"
            " --setRobustFitAlgo Minuit2"
            " --maxFailedSteps 100"
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
            f" --parallel {args.nparallel}"
            f" {args.combineargs}"
        )
        
        cmds.append(cmd)
        
        cmd = f"for dir in $(find {args.indir} -mindepth 0 -maxdepth 0 -type d); do plot1DScan.py $dir/{wkdir_scan}/higgsCombine.Test.MultiDimFit.mH120.root -o $dir/{wkdir_scan}/deltaNLL_r; done"
        cmds.append(cmd)
    
    if (args.impacts) :
        
        cmd = (
            "./run_impacts_limits.sh"
            f" \"{args.indir}\""
            f" {args.signal}"
        )
        
        cmds.append(cmd)
    
    if (args.gof) :
        
        #cmd = f"./run_gof_limits.sh \"{args.indir}\""
        #cmds.append(cmd)
        
        cmd = f"./run_gof_limits.sh \"{args.indir}\" b"
        cmds.append(cmd)
        
        cmd = f"./run_gof_limits.sh \"{args.indir}\" sb"
        cmds.append(cmd)
    
    #if (args.collgof) :
    #    
    #    cmd = (
    #        f"for dir in $(find {args.indir} -mindepth 0 -maxdepth 0 -type d); do \n"
    #        f"  nice -n 10 combineTool.py -M CollectGoodnessOfFit --input $dir/{wkdir_gof}/higgsCombine.Test.GoodnessOfFit.mH120.root $dir/{wkdir_gof}_b/higgsCombine.Test.GoodnessOfFit.root -o $dir/{wkdir_gof}_b/gof.json & \n"
    #        f"  nice -n 10 combineTool.py -M CollectGoodnessOfFit --input $dir/{wkdir_gof}/higgsCombine.Test.GoodnessOfFit.mH120.root $dir/{wkdir_gof}_sb/higgsCombine.Test.GoodnessOfFit.root -o $dir/{wkdir_gof}_sb/gof.json & \n"
    #        f"  while [[ $(jobs -rp | wc -l) -ge {args.nparallel} ]]; do \n"
    #        "       sleep 1 \n"
    #        "   done \n"
    #        " done \n"
    #        " wait $(jobs -rp) \n"
    #        " echo"
    #    )
    #    cmds.append(cmd)
    #    
    #    cmd = (
    #        f"for dir in $(find {args.indir} -mindepth 0 -maxdepth 0 -type d); do \n"
    #        f"  nice -n 10 plotGof.py $dir/{wkdir_gof}_b/gof.json --statistic saturated --mass 120.0 -o $dir/{wkdir_gof}_b/gof_plot_b --title-right=\"b-only\" & \n"
    #        f"  nice -n 10 plotGof.py $dir/{wkdir_gof}_sb/gof.json --statistic saturated --mass 120.0 -o $dir/{wkdir_gof}_sb/gof_plot_sb --title-right=\"s+b\" & \n"
    #        f"  while [[ $(jobs -rp | wc -l) -ge {args.nparallel} ]]; do \n"
    #        "       sleep 1 \n"
    #        "   done \n"
    #        " done \n"
    #        " wait $(jobs -rp) \n"
    #        " echo"
    #    )
    #    cmds.append(cmd)
    
    cmut.exec_cmds(cmds)
    
    
    if (args.colltoylimits) :
        
        d_mappings_toys = {
            "obs": ["obs", "", ""],
            "exp0": ["exp50", ".quant0.500", "--expectedFromGrid=0.5"],
            "exp-1": ["exp16", ".quant0.160", "--expectedFromGrid=0.16"],
            "exp+1": ["exp84", ".quant0.840", "--expectedFromGrid=0.84"],
            "exp-2": ["exp2p5", ".quant0.025", "--expectedFromGrid=0.025"],
            "exp+2": ["exp97p5", ".quant0.975", "--expectedFromGrid=0.975"],
        }
        
        d_error_root_files = {_key: [] for _key in d_mappings_toys.keys()}
        d_error_sub_files = {_key: [] for _key in d_mappings_toys.keys()}
        d_error_sub_files_missing = {_key: [] for _key in d_mappings_toys.keys()}
        
        for sig in cmut.natural_sort(glob.glob(f"{args.indir}*/")) :
            
            d_limits_tmp = {}
            
            for key, val in d_mappings_toys.items() :
                
                try :
                    
                    #fname_toylimits = f"{sig}/{wkdir_toylimits}_{val[0]}/condor/{args.timestamp}/higgsCombine.Test.HybridNew.mH120{val[1]}.root"
                    fname_toylimits = f"{sig}/{wkdir_toylimits}/condor/{args.timestamp}/higgsCombineTest.HybridNew.mH120{val[1]}.root"
                    rdf_toys = ROOT.RDataFrame("limit", fname_toylimits)
                    arr_limits_toys = rdf_toys.AsNumpy(["limit"])["limit"]
                    
                    d_limits_tmp[key] = arr_limits_toys[0]
                
                except Exception as excp :
                    
                    cmut.logger.warning(f"Error reading toy limits from [{fname_toylimits}]:\n{excp}")
                    d_error_root_files[key].append(fname_toylimits)
                    
                    fname_sub = os.path.join(os.path.dirname(fname_toylimits), f"condor_toylimits_{val[0]}.sub")
                    
                    if (os.path.isfile(fname_sub)) :
                        d_error_sub_files[key].append(fname_sub)
                    else :
                        d_error_sub_files_missing[key].append(fname_sub)
            d_limits_tmp = {"120.0": d_limits_tmp}
            outdir_tmp = f"{sig}/toylimits/{args.timestamp}"
            outfname_tmp = f"{outdir_tmp}/limits.json"
            
            os.system(f"mkdir -p {outdir_tmp}")
            
            with open(outfname_tmp, "w") as fopen:
                json.dump(d_limits_tmp, fopen, indent=4)
        
        n_error = sum([len(_val) for _val in d_error_root_files.values()])
        l_error_sub_files = list(itertools.chain(*d_error_sub_files.values()))
        l_error_sub_files_missing = list(itertools.chain(*d_error_sub_files_missing.values()))
        
        if (n_error) :
            
            print(f"Error reading toy limits from the following {n_error} files:")
            
            ifname = 0
            for key, flist in d_error_root_files.items() :
                
                for fname in flist :
                    
                    ifname += 1
                    n_chars = max([len(_key) for _key in d_error_root_files.keys()])
                    print(f"[{ifname:4d}] [{key:<{n_chars}}] {fname}")
            
            if (len(l_error_sub_files)) :
                fname = "error_toylimits_sub_files.txt"
                with open(fname, "w") as fopen :
                    fopen.write("\n".join(l_error_sub_files) + "\n")
                print("")
                cmut.logger.info(f"Written the list of {len(l_error_sub_files)} condor sub files to: {fname}")
            
            if (len(l_error_sub_files_missing)) :
                fname = "error_toylimits_sub_files_missing.txt"
                with open(fname, "w") as fopen :
                    fopen.write("\n".join(l_error_sub_files_missing) + "\n")
                print("")
                cmut.logger.info(f"Written the list of {len(l_error_sub_files_missing)} missing condor sub files to: {fname}")
    
    if l_toylimits_skipped :
        print(f"Skipped {len(l_toylimits_skipped)} toy limits:")
        for isig, (sig_key, indir, excp) in enumerate(l_toylimits_skipped) :
            print(f"[{isig+1:>4}] [{sig_key}] [{indir}]:\n{excp}")
    
    if (args.plotcorr) :
        
        l_fnames = glob.glob(f"{wkdir_plotcorr}/fitDiagnostics.Test.root")
        
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

