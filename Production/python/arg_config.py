import os

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing


############################## Parse arguments ##############################

d_procConfig = {
    "Data": {
        "2016_preVFP": {
            "condition": "auto:run2_data",
            "era": "Run2_2016_HIPM",
            "eramodifier": "run2_nanoAOD_106Xv2",
        },
        "2016_postVFP": {
            "condition": "auto:run2_data",
            "era": "Run2_2016",
            "eramodifier": "run2_nanoAOD_106Xv2",
        },
        "2017":{
            "condition": "auto:run2_data",
            "era": "Run2_2017",
            "eramodifier": "run2_nanoAOD_106Xv2",
        },
        "2018":{
            "condition": "auto:run2_data",
            "era": "Run2_2018",
            "eramodifier": "run2_nanoAOD_106Xv2",
        },
    },
    
    "MC": {
        "2016_preVFP": {
            "condition": "auto:run2_mc_pre_vfp",
            "era": "Run2_2016_HIPM",
            "eramodifier": "run2_nanoAOD_106Xv2",
        },
        "2016_postVFP": {
            "condition": "auto:run2_mc",
            "era": "Run2_2016",
            "eramodifier": "run2_nanoAOD_106Xv2",
        },
        "2017":{
            "condition": "auto:phase1_2017_realistic",
            "era": "Run2_2017",
            "eramodifier": "run2_nanoAOD_106Xv2",
        },
        "2018":{
            "condition": "auto:phase1_2018_realistic",
            "era": "Run2_2018",
            "eramodifier": "run2_nanoAOD_106Xv2",
        },
    }
}

d_procConfig["Embed"] = d_procConfig["Data"]

def get_args() :
    
    args = VarParsing("analysis")
    
    args.register("sourceFile",
        "", # Default value
        VarParsing.multiplicity.singleton, # singleton or list
        VarParsing.varType.string, # string, int, or float
        "File containing list of input files" # Description
    )

    # args.register("fileList",
    #     "[]",
    #     VarParsing.multiplicity.list,
    #     VarParsing.varType.string,
    #     "List of root files to process (alternative to sourceFile option).")

    args.register("outFile",
        "nanoaod.root", # Default value
        VarParsing.multiplicity.singleton, # singleton or list
        VarParsing.varType.string, # string, int, or float
        "Output file name" # Description
    )

    args.register("fileNamePrefix",
        "",
        VarParsing.multiplicity.singleton,
        VarParsing.varType.string,
        "Prefix to add to input file names.")

    args.register("lumiFile",
        "", 
        VarParsing.multiplicity.singleton, 
        VarParsing.varType.string,
        "JSON file with lumi mask.")
    
    args.register("eventRange",
        "", # Default value
        VarParsing.multiplicity.singleton, # singleton
        VarParsing.varType.string, # string, int, or float
        "Syntax: Run1:Event1-Run2:Event2 Run3:Event3-Run4:Event4(includes both)" # Description
    )
    
    args.register("sampleType",
        "", # Default value
        VarParsing.multiplicity.singleton, # singleton or list
        VarParsing.varType.string, # string, int, or float
        "Must be one of the following: [Data, MC, Embed]" # Description
    )
    
    args.register("era",
        "", # Default value
        VarParsing.multiplicity.singleton, # singleton or list
        VarParsing.varType.string, # string, int, or float
        "Must be one of the following: [2016, 2017, 2018]" # Description
    )
    
    args.register("disTauTagOutputOpt",
        0, # Default value
        VarParsing.multiplicity.singleton, # singleton or list
        VarParsing.varType.int, # string, int, or float
        "0: NANOAOD, 1: NANOAOD + disTauTagger, 2: disTauTagger only" # Description
    )
    
    args.register("pfCandExtraCut",
        "", # Default value
        VarParsing.multiplicity.singleton, # singleton or list
        VarParsing.varType.string, # string, int, or float
        "Enter NanoAOD compatible cut as a string" # Description
    )
    
    # args.register("debug",
    #     0, # Default value
    #     VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    #     VarParsing.VarParsing.varType.int, # string, int, or float
    #     "Print debug statements" # Description
    # )
    
    args.register("debugEDM",
        False, # Default value
        VarParsing.multiplicity.singleton, # singleton or list
        VarParsing.varType.bool, # string, int, or float
        "Create EDM file for debugging collection content" # Description
    )
    
    args.parseArguments()
    
    fNames = []
    
    if(len(args.inputFiles) and len(args.sourceFile)):
        raise ValueError("Error: fileList and sourceFile are interchangeable, only one should be specified.")

    if (len(args.inputFiles)) :
        fNames = args.inputFiles
    
    elif (len(args.sourceFile)) :
        with open(args.sourceFile) as f:
            args.inputFiles = f.readlines()
            args.inputFiles = [_fname.strip() for _fname in args.inputFiles]
            args.inputFiles = [_fname for _fname in args.inputFiles if _fname[0] != "#"]
    
    if ("/" in args.outFile) :
        outDir = args.outFile[0: args.outFile.rfind("/")]
        os.system("mkdir -p %s" %(outDir))
    
    if (args.sampleType not in d_procConfig.keys()) :
        raise ValueError(f"sampleType must be one of the following: [{', '.join(d_procConfig.keys())}]")
    
    if (args.era not in d_procConfig[args.sampleType].keys()) :
        raise ValueError(f"Era for {args.sampleType} must be one of the following: [{', '.join(d_procConfig[args.sampleType].keys())}]")
    
    return args
