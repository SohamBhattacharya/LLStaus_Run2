#!/usr/bin/env python3

import awkward
import coffea
import coffea.processor
import dataclasses
import glob
import hist
import logging
import numpy
import os
import sortedcontainers
import typing
import uproot

#import ROOT
#ROOT.gROOT.SetBatch(True)

#import CMS_lumi, tdrstyle
#import utils

from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

@dataclasses.dataclass
class MyProcessor(coffea.processor.ProcessorABC) :
    
    #datasets         : typing.List[str]
    dataset_args     : typing.Dict
    
    def __post_init__(self) :
        
        self.dataset_axis = hist.axis.StrCategory([], growth = True, name = "dataset", label = "dataset")
        
        self.d_hist_axis = {
            "pt"        : {"bins": 500, "start": 0, "stop": 5000}
        }
        
        self.d_hist_scheme = {
            "MET": {"pt": "pt"}
        }
        
        self._accumulator = {}
        
        for obj, d_qty in self.d_hist_scheme.items() :
            
            for qty, ax_key in d_qty.items() :
                
                d_ax_args = {}
                d_ax_args.update(self.d_hist_axis[ax_key])
                
                ax_name = f"{obj}_{qty}"
                d_ax_args["name"] = ax_name
                d_ax_args["label"] = ax_name
                
                self._accumulator[ax_name] = hist.Hist(self.dataset_axis, hist.axis.Regular(**d_ax_args))
        
        print(self._accumulator.keys())
    
    
    @property
    def accumulator(self) :
        
        return self._accumulator
    
    # we will receive a NanoEvents instead of a coffea DataFrame
    def process(self, events) :
        
        dataset = events.metadata["dataset"]
        args = self.dataset_args[dataset]
        
        filters = args["filters"]
        filter_str = " & ".join(filters)
        
        triggers = args["triggers"]
        trigger_str = " | ".join(triggers)
        
        events = events[eval(filter_str)]
        events = events[eval(trigger_str)]
        
        output = self.accumulator
        
        # Skip processing as it is an EmptyArray
        if not len(events) :
            
            return output
        
        for obj, d_qty in self.d_hist_scheme.items() :
            
            for qty, ax_key in d_qty.items() :
                
                d_ax_args = {}
                d_ax_args.update(self.d_hist_axis[ax_key])
                
                ax_name = f"{obj}_{qty}"
                
                d_hist_args = {"dataset": events.metadata["dataset"]}
                
                #if (qty not in events[obj]) :
                #    print(obj, events[obj].__dict__)
                    
                d_hist_args[ax_name] = awkward.flatten([events[obj][qty]], axis = None)
                output[ax_name].fill(**d_hist_args)
        
        return output
    
    
    def postprocess(self, accumulator):
        
        pass
        #return accumulator





def main() :
    
    #/MET/Run2017B-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD
    #/MET/Run2017C-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD
    #/MET/Run2017D-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD
    #/MET/Run2017E-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD
    #/MET/Run2017F-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD
    
    #/pnfs/desy.de/cms/tier2/store/data/Run2017B/MET/NANOAOD/UL2017_MiniAODv2_NanoAODv9-v1/*/*.root
    #/pnfs/desy.de/cms/tier2/store/data/Run2017D/MET/NANOAOD/UL2017_MiniAODv2_NanoAODv9-v1/*/*.root
    #/pnfs/desy.de/cms/tier2/store/data/Run2017E/MET/NANOAOD/UL2017_MiniAODv2_NanoAODv9-v1/*/*.root
    #/pnfs/desy.de/cms/tier2/store/data/Run2017F/MET/NANOAOD/UL2017_MiniAODv2_NanoAODv9-v1/*/*.root

    base_storage_dir = "/pnfs/desy.de/cms/tier2/store/user/sobhatta/LongLivedStaus/NanoAOD"
    
    datasets = {
        "MET_Run2017B_UL":      f"{base_storage_dir}/2017/MET/crab_MET_Run2017B_UL/240710_012134/*/*.root",
        "MET_Run2017C_UL":      f"{base_storage_dir}/2017/MET/crab_MET_Run2017C_UL/240710_012141/*/*.root",
        "MET_Run2017D_UL":      f"{base_storage_dir}/2017/MET/crab_MET_Run2017D_UL/240710_085120/*/*.root",
        "MET_Run2017E_UL":      f"{base_storage_dir}/2017/MET/crab_MET_Run2017E_UL/240710_012147/*/*.root",
        "MET_Run2017F_UL":      f"{base_storage_dir}/2017/MET/crab_MET_Run2017F_UL/240710_012153/*/*.root",
    }
    
    d_filters = {}
    d_filters["2017UL"] = [
        "events.Flag.goodVertices",
        "events.Flag.globalSuperTightHalo2016Filter",
        "events.Flag.HBHENoiseFilter",
        "events.Flag.HBHENoiseIsoFilter",
        "events.Flag.EcalDeadCellTriggerPrimitiveFilter",
        "events.Flag.BadPFMuonFilter",
        "events.Flag.BadPFMuonDzFilter",
        "events.Flag.eeBadScFilter",
        "events.Flag.ecalBadCalibFilter",
    ]
    
    d_triggers = {}
    d_triggers["2017UL"] = {
        "events.HLT.PFMET120_PFMHT120_IDTight",
        "events.HLT.PFMET130_PFMHT130_IDTight",
        "events.HLT.PFMET140_PFMHT140_IDTight"
    }
    
    dataset_args = {
        "MET_Run2017B_UL":      {"era": "2017UL"},
        "MET_Run2017C_UL":      {"era": "2017UL"},
        "MET_Run2017D_UL":      {"era": "2017UL"},
        "MET_Run2017E_UL":      {"era": "2017UL"},
        "MET_Run2017F_UL":      {"era": "2017UL"},
    }
    
    for key, val in dataset_args.items() :
        
        era = val["era"]
        val["filters"] = d_filters[era]
        val["triggers"] = d_triggers[era]
    
    for key, val in datasets.items() :
        
        l_file_tmp = glob.glob(val)#[0: 20]
        datasets[key] = [f"{_ele}" for _ele in l_file_tmp]
    
    #print(datasets)
    
    output = coffea.processor.run_uproot_job(
        datasets,
        "Events",
        MyProcessor(
            dataset_args = dataset_args,
        ),
        #executor = coffea.processor.iterative_executor,
        executor = coffea.processor.futures_executor,
        executor_args = {
            "schema": NanoAODSchema,
            "skipbadfiles": True,
            #"xrootdtimeout": 600, #sec
            "workers": 15,
        },
        #chunksize=1000
    )
    
    print(output)
    outfile = "output/met_study/met_study.root"
    outdir = os.path.dirname(outfile)
    os.system(f"mkdir -p {outdir}")
    
    with uproot.recreate(outfile) as fout:
        
        for dataset_key in datasets.keys() :
            
            for hist_key, histo in output.items() :
                
                for ax in histo.axes :
                    
                    if ax.name == "dataset" :
                        continue
                    
                    print(dataset_key, hist_key, ax)
                    h_tmp = histo[{"dataset": dataset_key}].project(ax.name)
                    h_tmp = h_tmp / h_tmp.sum(flow = True)
                    fout[f"{dataset_key}/{ax.name}"] = h_tmp
    
    return 0


if __name__ == "__main__" :
    
    main()
