#!/usr/bin/env python3

import argparse
import sortedcontainers
import ROOT

import utils.commonutils as cmut

def main() :
    
    # Argument parser
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        "--dir",
        help = " ".join([
            "Input directory containing SMS-TStauStau_MStau-*.",
            "For e.g. results/limits_signal_v12/llstau_maximally-mixed/channels_BRT2/eras_all/",
        ]),
        type = str,
        required = True,
    )
    
    parser.add_argument(
        "--toystamp",
        help = "Timestamp (usually YYYY-MM-DD_hh-mm-ss) of the toy condor jobs",
        type = str,
        required = True,
    )
    
    parser.add_argument(
        "--limits",
        help = "Limits to compare",
        type = str,
        nargs = "+",
        required = False,
        default = [
            "obs",
            "exp0",
            "exp-1",
            "exp+1",
            "exp-2",
            "exp+2",
        ]
    )
    
    
    # Parse arguments
    args = parser.parse_args()
    
    l_sigs = [
        #"SMS-TStauStau_MStau-100_ctau-3mm_mLSP-1_from_ctau-5mm",
        #"SMS-TStauStau_MStau-100_ctau-50mm_mLSP-1",
        #"SMS-TStauStau_MStau-100_ctau-100mm_mLSP-1",
        #"SMS-TStauStau_MStau-100_ctau-300mm_mLSP-1_from_ctau-1000mm",
        #
        #"SMS-TStauStau_MStau-200_ctau-3mm_mLSP-1_from_ctau-5mm",
        #"SMS-TStauStau_MStau-200_ctau-50mm_mLSP-1",
        #"SMS-TStauStau_MStau-200_ctau-100mm_mLSP-1",
        #"SMS-TStauStau_MStau-200_ctau-300mm_mLSP-1_from_ctau-1000mm",
        #
        #"SMS-TStauStau_MStau-300_ctau-3mm_mLSP-1_from_ctau-5mm",
        #"SMS-TStauStau_MStau-300_ctau-50mm_mLSP-1",
        #"SMS-TStauStau_MStau-300_ctau-100mm_mLSP-1",
        #"SMS-TStauStau_MStau-300_ctau-300mm_mLSP-1_from_ctau-1000mm",
        #
        #"SMS-TStauStau_MStau-400_ctau-3mm_mLSP-1_from_ctau-5mm",
        #"SMS-TStauStau_MStau-400_ctau-50mm_mLSP-1",
        #"SMS-TStauStau_MStau-400_ctau-100mm_mLSP-1",
        #"SMS-TStauStau_MStau-400_ctau-300mm_mLSP-1_from_ctau-1000mm",
        
        
        # ATLAS points
        "SMS-TStauStau_MStau-100_ctau-3mm_mLSP-1_from_ctau-5mm",
        "SMS-TStauStau_MStau-100_ctau-30mm_mLSP-1_from_ctau-50mm",
        "SMS-TStauStau_MStau-100_ctau-300mm_mLSP-1_from_ctau-1000mm",
        
        "SMS-TStauStau_MStau-200_ctau-3mm_mLSP-1_from_ctau-5mm",
        "SMS-TStauStau_MStau-200_ctau-30mm_mLSP-1_from_ctau-50mm",
        "SMS-TStauStau_MStau-200_ctau-300mm_mLSP-1_from_ctau-1000mm",
        
        "SMS-TStauStau_MStau-300_ctau-3mm_mLSP-1_from_ctau-5mm",
        "SMS-TStauStau_MStau-300_ctau-30mm_mLSP-1_from_ctau-50mm",
        "SMS-TStauStau_MStau-300_ctau-300mm_mLSP-1_from_ctau-1000mm",
        
        "SMS-TStauStau_MStau-400_ctau-3mm_mLSP-1_from_ctau-5mm",
        "SMS-TStauStau_MStau-400_ctau-30mm_mLSP-1_from_ctau-50mm",
        "SMS-TStauStau_MStau-400_ctau-300mm_mLSP-1_from_ctau-1000mm",
    ]
    
    d_limits = {}
    #d_limits = sortedcontainers.SortedDict()
    
    d_mappings_asymptotic = {
        "obs": 5,
        "exp0": 2,
        "exp-1": 1,
        "exp+1": 3,
        "exp-2": 0,
        "exp+2": 4,
    }
    
    d_mappings_toys = {
        "obs": ["obs", ""],
        "exp0": ["exp50", ".quant0.500"],
        "exp-1": ["exp16", ".quant0.160"],
        "exp+1": ["exp84", ".quant0.840"],
        "exp-2": ["exp2p5", ".quant0.025"],
        "exp+2": ["exp97p5", ".quant0.975"],
    }
    
    for sig in l_sigs :
        
        sig_info = cmut.parse_string_regex(
            s = sig,
            regexp = "SMS-TStauStau_MStau-(?P<mstau>\w+)_ctau-(?P<ctau>\w+)mm_mLSP-(?P<mlsp>\d+)",
        )
        
        sig_key = (sig_info["mstau"], sig_info["ctau"])
        
        if sig_key not in d_limits :
            
            d_limits[sig_key] = {
                "asymptotic": {_lim: 0 for _lim in args.limits},
                "toys": {_lim: 0 for _lim in args.limits},
            }
        
        rdf_asymptotic = ROOT.RDataFrame("limit", f"{args.dir}/{sig}/limits/higgsCombine.Test.AsymptoticLimits.mH120.root")
        arr_limits_asymptotic = rdf_asymptotic.AsNumpy(["limit"])["limit"]
        
        for key, val in d_mappings_asymptotic.items() :
            
            if key not in d_limits[sig_key]["asymptotic"] :
                continue
            
            d_limits[sig_key]["asymptotic"][key] = arr_limits_asymptotic[val]
        
        for key, val in d_mappings_toys.items() :
            
            if key not in d_limits[sig_key]["toys"] :
                continue
            
            #rdf_toys = ROOT.RDataFrame("limit", f"{args.dir}/{sig}/toylimits_{val[0]}/higgsCombine.Test.HybridNew.mH120{val[1]}.root")
            #rdf_toys = ROOT.RDataFrame("limit", f"{args.dir}/{sig}/toylimits_{val[0]}/condor/{args.toystamp}/higgsCombine.Test.HybridNew.mH120{val[1]}.root")
            rdf_toys = ROOT.RDataFrame("limit", f"{args.dir}/{sig}/toylimits/condor/{args.toystamp}/higgsCombineTest.HybridNew.mH120{val[1]}.root")
            arr_limits_toys = rdf_toys.AsNumpy(["limit"])["limit"]
            
            d_limits[sig_key]["toys"][key] = arr_limits_toys[0]
    
    # Print results
    print()
    
    header = [
        "mass [GeV]",
        "ctau0 [mm]",
        "limit",
        "asymp.",
        "toys",
        #"asymp./toys",
        "toys/asymp.",
    ]
    
    #str_fmt = "".join([f"{{:<{max(6, len(_x))}}}" for _x in header])
    #str_fmt = "".join([f"{{:<{max(6, len(_x)+2)}}}" for _x in header])
    str_fmt = ",".join([f"{{:<{len(_x)+2}}}" for _x in header])
    
    print(str_fmt.format(*header))
    print()
    
    for sig_key, limits in d_limits.items() :
        
        limits_asymptotic = limits["asymptotic"]
        limits_toys = limits["toys"]
        
        for lim_key in limits_asymptotic.keys() :
            
            row = [
                f"{sig_key[0]}",
                f"{sig_key[1]}",
                f"{lim_key}",
                f"{limits_asymptotic[lim_key]:0.2f}",
                f"{limits_toys[lim_key]:0.2f}",
                #f"{limits_asymptotic[lim_key]/limits_toys[lim_key]:0.2f}",
                f"{limits_toys[lim_key]/limits_asymptotic[lim_key]:0.2f}",
            ]
            
            print(str_fmt.format(*row))
        
        print()
    
    return 0

if __name__ == "__main__":
    main()