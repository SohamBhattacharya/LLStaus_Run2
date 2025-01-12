#!/usr/bin/env python3

import awkward
import coffea
import coffea.processor
import dataclasses
import hist
import os
import sortedcontainers
import typing
import uproot

import ROOT
ROOT.gROOT.SetBatch(True)

from coffea.nanoevents import NanoEventsFactory, NanoAODSchema


@dataclasses.dataclass
class MyProcessor(coffea.processor.ProcessorABC) :
    
    #datasets         : typing.List[str]
    dataset_args     : typing.Dict
    
    def __post_init__(self) :
        
        #self.dataset_axis = coffea.hist.Cat("dataset", "dataset")
        
        #self.dataset_axis = hist.axis.StrCategory(self.datasets, growth = True, name = "dataset", label = "dataset")
        #
        #self._accumulator = sortedcontainers.SortedDict({
        #    "GenVisTauh": hist.Hist(
        #        self.dataset_axis,
        #        hist.axis.Regular(100, 0, 1000, name = "GenVisTauh_e", label = "GenVisTauh_e"),
        #        hist.axis.Regular(100, 0, 1000, name = "GenVisTauh_pt", label = "GenVisTauh_pt"),
        #        hist.axis.Regular(200, 0, 2, name = "GenVisTauh_e_by_stau_mass", label = "GenVisTauh_e_by_stau_mass"),
        #        storage = "weight",
        #        name = "Counts"
        #    ),
        #    
        #    "GenVisTaul": hist.Hist(
        #        self.dataset_axis,
        #        hist.axis.Regular(100, 0, 1000, name = "GenVisTaul_e", label = "GenVisTaul_e"),
        #        hist.axis.Regular(100, 0, 1000, name = "GenVisTaul_pt", label = "GenVisTaul_pt"),
        #        hist.axis.Regular(200, 0, 2, name = "GenVisTaul_e_by_stau_mass", label = "GenVisTaul_e_by_stau_mass"),
        #        storage = "weight",
        #        name = "Counts"
        #    ),
        #})
        
        self.dataset_axis = hist.axis.StrCategory([], growth = True, name = "dataset", label = "dataset")
        
        self.d_hist_axis = {
            "e": {"bins": 100, "start": 0, "stop": 1000},
            "pt": {"bins": 100, "start": 0, "stop": 1000},
            "ebym": {"bins": 200, "start": 0, "stop": 2},
        }
        
        self.d_hist_scheme = {
            "GenVisTauh": {"pt": "pt", "ebym": "ebym"},
            #"GenVisTauh": {"energy": "e", "pt": "pt", "ebym": "ebym"},
            #"GenVisTaul": {"energy": "e", "pt": "pt", "ebym": "ebym"},
            
            "GenVisTauh1": {"pt": "pt"},
            "GenVisTauh2": {"pt": "pt"},
            
            "Jet1": {"pt": "pt"},
            "Jet2": {"pt": "pt"},
            
            "MET": {"pt": "pt"},
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
        
        output = self.accumulator
        
        stau_mass = 250.0
        
        sel_idx = (awkward.num(events.GenVisTau, axis = 1) == 2)
        events = events[sel_idx]
        
        GenStau = events.GenPart[
            ((abs(events.GenPart.pdgId) == 1000015) | (abs(events.GenPart.pdgId) == 2000015))
            #& (abs(events.GenPart.pdgId[events.GenPart.genPartIdxMother]) != 1000015)
            & events.GenPart.hasFlags(["isHardProcess"])
            & events.GenPart.hasFlags(["isFirstCopy"])
            & (events.GenPart.mass == stau_mass)
        ]
        
        GenLsp = events.GenPart[
            (abs(events.GenPart.pdgId) == 1000022)
            #& (abs(events.GenPart.pdgId[events.GenPart.genPartIdxMother]) == 1000015)
            & events.GenPart.hasFlags(["isHardProcess"])
            & events.GenPart.hasFlags(["isFirstCopy"])
            & (events.GenPart.mass == 1)
        ]
        
        sel_idx = (awkward.num(GenStau, axis = 1) == 2) & (awkward.num(GenLsp, axis = 1) == 2)
        events = events[sel_idx]
        
        #GenTau = events.GenPart[
        #    (abs(events.GenPart.pdgId) == 15)
        #    & events.GenPart.hasFlags(["isHardProcess"])
        #    & events.GenPart.hasFlags(["isFirstCopy"])
        #]
        
        
        #GenVisTaul = events.GenPart[
        #    ((abs(events.GenPart.pdgId) == 11) | (abs(events.GenPart.pdgId) == 13))
        #    & events.GenPart.hasFlags(["isFirstCopy"])
        #    & events.GenPart.hasFlags(["isDirectHardProcessTauDecayProduct"])
        #]
        #
        #print(GenVisTaul)
        #print(awkward.num(GenVisTaul, axis = 1))
        #
        #GenStaul = GenVisTaul.distinctParent.distinctParent
        #GenVisTaul_stauRF = GenVisTaul.boost(-GenStaul.boostvec)
        
        
        events["GenVisTauh"] = events.GenVisTau
        
        GenVisTauh_sorted = events.GenVisTauh[awkward.argsort(events.GenVisTauh.pt, axis = 1, ascending = False)]
        events["GenVisTauh1"] = GenVisTauh_sorted[:, 0]
        events["GenVisTauh2"] = GenVisTauh_sorted[:, 1]
        
        events["GenTauh"] = events.GenVisTau.parent
        
        #GenStauh = events.GenTauh.distinctParent
        #GenStauh = GenStauh[
        #    ((abs(GenStauh.pdgId) == 1000015) | (abs(GenStauh.pdgId) == 2000015))
        #    | ((abs(GenStauh.distinctParent.pdgId) == 1000015) | (abs(GenStauh.distinctParent.pdgId) == 2000015))
        #]
        
        GenStauh_set0 = events.GenTauh
        GenStauh_set0 = GenStauh_set0[((abs(GenStauh_set0.pdgId) == 1000015) | (abs(GenStauh_set0.pdgId) == 2000015))]
        #GenStauh_set0[~awkward.is_none(GenStauh_set0)]
        
        GenStauh_set1 = events.GenTauh.distinctParent
        GenStauh_set1 = GenStauh_set1[((abs(GenStauh_set1.pdgId) == 1000015) | (abs(GenStauh_set1.pdgId) == 2000015))]
        #GenStauh_set1[~awkward.is_none(GenStauh_set1)]
        
        GenStauh_set2 = events.GenTauh.distinctParent.distinctParent
        GenStauh_set2 = GenStauh_set2[((abs(GenStauh_set2.pdgId) == 1000015) | (abs(GenStauh_set2.pdgId) == 2000015))]
        
        #GenStauh_set2[~awkward.is_none(GenStauh_set2)]
        GenStauh = awkward.concatenate([GenStauh_set0, GenStauh_set1, GenStauh_set2], axis = 1)
        GenStauh = GenStauh[~awkward.is_none(GenStauh, axis = 1)]
        
        #events["GenStauh"] = GenStauh_set1
        events["GenStauh"] = GenStauh
        
        #print(events.GenVisTauh)
        #print(events.GenStauh)
        #print(awkward.num(events.GenStauh, axis = 1))
        #print(events.GenStauh.boostvec)
        
        sel_idx = (awkward.num(events.GenStauh, axis = 1) == 2)
        events = events[sel_idx]
        
        events["GenVisTauh_stauRF"] = events.GenVisTauh.boost(-events.GenStauh.boostvec)
        events["GenVisTauh", "ebym"] = events.GenVisTauh_stauRF.energy / events.GenStauh.mass
        
        #print("GenTauh.pdgId:", events.GenTauh.pdgId)
        #print("GenStauh.pdgId:", events.GenStauh.pdgId)
        ##print("GenStauh.distinctParent.pdgId:", GenStauh.distinctParent.pdgId)
        #print("num(events.GenVisTau):", awkward.num(events.GenVisTau, axis = 1))
        #print("GenVisTauh.ebym:", events.GenVisTauh.ebym)
        
        jets = events.Jet[
            (events.Jet.pt > 30)
            & (abs(events.Jet.eta) < 2.4)
            & (events.Jet.jetId >= 6)
        ]
        
        (dr, (_, _jets)) = jets.metric_table(events.GenVisTauh, metric = coffea.nanoevents.methods.vector.LorentzVector.delta_r, return_combinations = True)
        jets = _jets[(dr < 0.4)]
        jets = jets[~awkward.is_none(jets, axis = 1)]
        
        events["Jet"] = jets
        sel_idx = (awkward.num(events.Jet, axis = 1) >= 2)
        events = events[sel_idx]
        
        events["Jet1"] = events.Jet[:, 0]
        events["Jet2"] = events.Jet[:, 1]
        
        # Skip processing as it is an EmptyArray
        if not len(events) :
            
            return output
        
        #output["GenVisTauh"].fill(
        #    dataset = events.metadata["dataset"],
        #    GenVisTauh_e = awkward.flatten(events.GenVisTau.energy),
        #    GenVisTauh_pt = awkward.flatten(events.GenVisTau.pt),
        #    GenVisTauh_e_by_stau_mass = awkward.flatten(events.GenVisTau.pt),
        #)
        
        #output["GenVisTaul"].fill(
        #    dataset = events.metadata["dataset"],
        #    GenVisTaul_e = awkward.flatten(GenVisTaul.energy),
        #    GenVisTaul_pt = awkward.flatten(GenVisTaul.pt),
        #    GenVisTaul_e_by_stau_mass = awkward.flatten(GenVisTaul_stauRF.energy / GenStaul.mass),
        #)
        
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
        
        return accumulator





def main() :
    
    d_fnamelist = {}
    #d_fnamelist["stau_LH"] = list(numpy.loadtxt("../Production/configs/sourceFiles/SMS-TStauStau_lefthanded_mStau-225to250_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv7-Nano02Apr2020_GridpackScan_102X_upgrade2018_realistic_v21-v1_NANOAODSIM/SMS-TStauStau_lefthanded_mStau-225to250_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv7-Nano02Apr2020_GridpackScan_102X_upgrade2018_realistic_v21-v1_NANOAODSIM.txt", dtype = str))
    #d_fnamelist["stau_RH"] = list(numpy.loadtxt("../Production/configs/sourceFiles/SMS-TStauStau_righthanded_mStau-225to250_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv7-Nano02Apr2020_GridpackScan_102X_upgrade2018_realistic_v21-v1_NANOAODSIM/SMS-TStauStau_righthanded_mStau-225to250_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv7-Nano02Apr2020_GridpackScan_102X_upgrade2018_realistic_v21-v1_NANOAODSIM.txt", dtype = str))
    #d_fnamelist["stau_MM"] = list(numpy.loadtxt("../Production/configs/sourceFiles/SMS-TStauStau_ctau-0p01to10_mStau-250to500_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv7-Nano02Apr2020_GridpackScan_102X_upgrade2018_realistic_v21-v1_NANOAODSIM/SMS-TStauStau_ctau-0p01to10_mStau-250to500_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv7-Nano02Apr2020_GridpackScan_102X_upgrade2018_realistic_v21-v1_NANOAODSIM.txt", dtype = str))
    
    #d_fnamelist["stau_LH"] = ["/pnfs/desy.de/cms/tier2/store/user/sobhatta/LongLivedStaus/NanoAOD/SUS-RunIISummer20UL18GEN-stau250_lsp1_ctau100mm_v6/crab_stau250_lsp1_ctau100mm/230311_014339/0000/nanoaod_with-disTauTagScore_1.root"]
    
    d_fnamelist["stau_LH"] = ["SMS-TStauStau_lefthanded_mStau-225to250_TuneCP5_13TeV-madgraphMLM-pythia8.root"]
    d_fnamelist["stau_RH"] = ["SMS-TStauStau_righthanded_mStau-225to250_TuneCP5_13TeV-madgraphMLM-pythia8.root"]
    d_fnamelist["stau_MM"] = ["SMS-TStauStau_ctau-0p01to10_mStau-250to500_TuneCP5_13TeV-madgraphMLM-pythia8.root"]
    #print(d_fnamelist["stau_LH"].shape)
    #print(d_fnamelist["stau_LH"][0: 1])
    
    datasets = sortedcontainers.SortedDict({
        "stau_LH": d_fnamelist["stau_LH"],
        "stau_RH": d_fnamelist["stau_RH"],
        "stau_MM": d_fnamelist["stau_MM"],
    })
    
    dataset_args = {
        "stau_LH": {},
        "stau_RH": {},
        "stau_MM": {},
    }
    
    #output = coffea.processor.run_uproot_job(
    #    datasets,
    #    "Events",
    #    MyProcessor(
    #        datasets = list(datasets.keys())
    #    ),
    #    #executor = coffea.processor.iterative_executor,
    #    executor = coffea.processor.futures_executor,
    #    executor_args = {
    #        "schema": NanoAODSchema,
    #        #"skipbadfiles": True,
    #        "xrootdtimeout": 600, #sec
    #        "workers": 10
    #    },
    #)
    
    output = coffea.processor.run_uproot_job(
        datasets,
        "Events",
        MyProcessor(
            dataset_args = dataset_args,
        ),
        executor = coffea.processor.futures_executor,
        executor_args = {
            "schema": NanoAODSchema,
            "skipbadfiles": True,
            "workers": 10,
        },
    )
    
    print(output)
    
    #with uproot.recreate("output_stau_kinematics_study.root") as fout:
    #    
    #    for key in output:
    #        
    #        histo = output[key]
    #        print(key)
    #        print(histo.axes)
    #        #print(histo.axes["dataset"])
    #        print(histo.__dict__)
    #        print(histo.axes.__dict__)
    #        
    #        for s in histo.axes["dataset"]:
    #            
    #            for ax in histo.axes :
    #                
    #                if ax.name == "dataset" :
    #                    
    #                    continue
    #                    
    #                print(s, ax)
    #                fout[f"{s}/{ax.name}"] = histo[{"dataset": s}].project(ax.name)
    
    
    with uproot.recreate("output/stau_kinematics_study/output_stau_kinematics_study.root") as fout:
        
        for dataset_key in datasets.keys() :
            
            for hist_key, histo in output.items() :
                
                for ax in histo.axes :
                    
                    if ax.name == "dataset" :
                        continue
                    
                    print(dataset_key, ax)
                    h_tmp = histo[{"dataset": dataset_key}].project(ax.name)
                    h_tmp = h_tmp / h_tmp.sum(flow = True)
                    fout[f"{dataset_key}/{ax.name}"] = h_tmp
    
    return 0


if __name__ == "__main__" :
    
    main()



#/SMS-TStauStau_lefthanded_mStau-225to250_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_GridpackScan_102X_upgrade2018_realistic_v21-v1/NANOAODSIM
#/SMS-TStauStau_righthanded_mStau-225to250_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_GridpackScan_102X_upgrade2018_realistic_v21-v1/NANOAODSIM
#/SMS-TStauStau_ctau-0p01to10_mStau-250to500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_GridpackScan_102X_upgrade2018_realistic_v21-v1/NANOAODSIM

#/SMS-TStauStau_ctau-0p01to10_mStau-90_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_GridpackScan_102X_upgrade2018_realistic_v21-v1/NANOAODSIM