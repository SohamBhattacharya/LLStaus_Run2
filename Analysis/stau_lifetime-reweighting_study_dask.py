#!/usr/bin/env python3

import awkward
import dask_awkward
import coffea
import coffea.processor
import coffea.dataset_tools
import dask
import dask_awkward
import dataclasses
import glob
import hist
import hist.dask
import logging
import numpy
import sortedcontainers
import typing
import uproot

#import ROOT
#ROOT.gROOT.SetBatch(True)

#import CMS_lumi, tdrstyle
#import utils

from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from dask.diagnostics import ProgressBar


@dataclasses.dataclass
class MyProcessor(coffea.processor.ProcessorABC) :
    
    #datasets         : typing.List[str]
    dataset_args     : typing.Dict
    
    def __post_init__(self) :
        
        #self.dataset_axis = hist.axis.StrCategory(self.datasets, growth = True, name = "dataset", label = "dataset")
        self.dataset_axis = hist.axis.StrCategory([], growth = True, name = "dataset", label = "dataset")
        
        #axis_pt = hist.axis.Regular(100, 0, 1000, name = "GenStau_pt", label = "GenStau_pt")
        #axis_path = hist.axis.Regular(1000, 0, 100, name = "GenStau_pathL", label = "GenStau_pathL")
        #axis_pxy = hist.axis.Regular(200, -1000, 1000, name = "MET_px", label = "MET_px"),
        
        self.d_hist_axis = {
            "pt"        : {"bins": 100, "start": 0, "stop": 1000},
            "pxy"       : {"bins": 200, "start": -1000, "stop": 1000},
            "pathL"     : {"bins": 1000, "start": 0, "stop": 100},
            "score"     : {"bins": 100, "start": 0, "stop": 1},
            "dxy"       : {"bins": 2000, "start": -100, "stop": 100},
            
            #"pathL"      : hist.axis.Regular(1000, 0, 100, name = "pathL", label = "pathL"),
            #"cpathL"      : hist.axis.Regular(1000, 0, 100, name = "pathL", label = "pathL"),
            #"px"        : hist.axis.Regular(200, -1000, 1000, name = "px", label = "px"),
            #"pt"        : hist.axis.Regular(200, -1000, 1000, name = "py", label = "py"),
        }
        
        self.d_hist_scheme = {
            #"GenStau": {"pt": "pt", "pathL": "pathL", "cpathL": "pathL"},
            "GenStau1": {"pt": "pt", "pathL": "pathL", "cpathL": "pathL"},
            "GenStau2": {"pt": "pt", "pathL": "pathL", "cpathL": "pathL"},
            
            #"GenTau": {"pt": "pt", "pathL": "pathL", "cpathL": "pathL", "vertexR": "pathL"},
            "GenTau1": {"pt": "pt", "vertexR": "pathL"},
            "GenTau2": {"pt": "pt", "vertexR": "pathL"},
            
            "Jet1": {"pt": "pt", "disTauTag_score1": "score", "dxy": "dxy"},
            "Jet2": {"pt": "pt", "disTauTag_score1": "score", "dxy": "dxy"},
            
            "MET": {"pt": "pt"}#, "px": "pxy", "py": "pxy"},
        }
        
        self._accumulator = {}
        
        for obj, d_qty in self.d_hist_scheme.items() :
            
            for qty, ax_key in d_qty.items() :
                
                d_ax_args = {}
                d_ax_args.update(self.d_hist_axis[ax_key])
                
                ax_name = f"{obj}_{qty}"
                d_ax_args["name"] = ax_name
                d_ax_args["label"] = ax_name
                
                self._accumulator[ax_name] = hist.dask.Hist(self.dataset_axis, hist.axis.Regular(**d_ax_args))
        
        #self._accumulator = {
        #    #"GenStau": hist.dask.Hist(
        #    #    self.dataset_axis,
        #    #    hist.axis.Regular(100, 0, 1000, name = "GenStau_pt", label = "GenStau_pt"),
        #    #    hist.axis.Regular(1000, 0, 100, name = "GenStau_pathL", label = "GenStau_pathL"),
        #    #    hist.axis.Regular(1000, 0, 100, name = "GenStau_cpathL", label = "GenStau_cpathL"),
        #    #    storage = "weight",
        #    #    name = "Counts"
        #    #),
        #    
        #    "GenStau_pt": hist.dask.Hist(self.dataset_axis, hist.axis.Regular(100, 0, 1000, name = "GenStau_pt", label = "GenStau_pt")),
        #    "GenStau_pathL": hist.dask.Hist(self.dataset_axis, hist.axis.Regular(1000, 0, 100, name = "GenStau_pathL", label = "GenStau_pathL")),
        #    "GenStau_cpathL": hist.dask.Hist(self.dataset_axis, hist.axis.Regular(1000, 0, 100, name = "GenStau_cpathL", label = "GenStau_cpathL")),
        #    
        #    
        #    #"GenStau1": hist.dask.Hist(
        #    #    self.dataset_axis,
        #    #    hist.axis.Regular(100, 0, 1000, name = "GenStau1_pt", label = "GenStau1_pt"),
        #    #    hist.axis.Regular(1000, 0, 100, name = "GenStau1_pathL", label = "GenStau1_pathL"),
        #    #    hist.axis.Regular(1000, 0, 100, name = "GenStau1_cpathL", label = "GenStau1_cpathL"),
        #    #    storage = "weight",
        #    #    name = "Counts"
        #    #),
        #    #
        #    #"GenStau2": hist.dask.Hist(
        #    #    self.dataset_axis,
        #    #    hist.axis.Regular(100, 0, 1000, name = "GenStau2_pt", label = "GenStau2_pt"),
        #    #    hist.axis.Regular(1000, 0, 100, name = "GenStau2_pathL", label = "GenStau2_pathL"),
        #    #    hist.axis.Regular(1000, 0, 100, name = "GenStau2_cpathL", label = "GenStau2_cpathL"),
        #    #    storage = "weight",
        #    #    name = "Counts"
        #    #),
        #    #
        #    #"GenTau": hist.dask.Hist(
        #    #    self.dataset_axis,
        #    #    hist.axis.Regular(100, 0, 1000, name = "GenTau_pt", label = "GenTau_pt"),
        #    #    hist.axis.Regular(1000, 0, 100, name = "GenTau_vertexR", label = "GenTau_vertexR"),
        #    #    storage = "weight",
        #    #    name = "Counts"
        #    #),
        #    #
        #    #"GenTau1": hist.dask.Hist(
        #    #    self.dataset_axis,
        #    #    hist.axis.Regular(100, 0, 1000, name = "GenTau1_pt", label = "GenTau1_pt"),
        #    #    hist.axis.Regular(1000, 0, 100, name = "GenTau1_vertexR", label = "GenTau1_vertexR"),
        #    #    storage = "weight",
        #    #    name = "Counts"
        #    #),
        #    #
        #    #"GenTau2": hist.dask.Hist(
        #    #    self.dataset_axis,
        #    #    hist.axis.Regular(100, 0, 1000, name = "GenTau2_pt", label = "GenTau2_pt"),
        #    #    hist.axis.Regular(1000, 0, 100, name = "GenTau2_vertexR", label = "GenTau2_vertexR"),
        #    #    storage = "weight",
        #    #    name = "Counts"
        #    #),
        #    #
        #    #"Jet1": hist.dask.Hist(
        #    #    self.dataset_axis,
        #    #    hist.axis.Regular(100, 0, 1000, name = "Jet1_pt", label = "Jet1_pt"),
        #    #    hist.axis.Regular(100, 0, 1, name = "Jet1_disTauTag_score1", label = "Jet1_disTauTag_score1"),
        #    #    storage = "weight",
        #    #    name = "Counts"
        #    #),
        #    #
        #    #"Jet2": hist.dask.Hist(
        #    #    self.dataset_axis,
        #    #    hist.axis.Regular(100, 0, 1000, name = "Jet2_pt", label = "Jet2_pt"),
        #    #    hist.axis.Regular(100, 0, 1, name = "Jet2_disTauTag_score1", label = "Jet2_disTauTag_score1"),
        #    #    storage = "weight",
        #    #    name = "Counts"
        #    #),
        #    #
        #    #"MET": hist.dask.Hist(
        #    #    self.dataset_axis,
        #    #    hist.axis.Regular(100, 0, 1000, name = "MET_pt", label = "MET_pt"),
        #    #    hist.axis.Regular(200, -1000, 1000, name = "MET_px", label = "MET_px"),
        #    #    hist.axis.Regular(200, -1000, 1000, name = "MET_py", label = "MET_py"),
        #    #    storage = "weight",
        #    #    name = "Counts"
        #    #),
        #}
        
        print(self._accumulator.keys())
    
    
    @property
    def accumulator(self) :
        
        return self._accumulator
    
    def pfcand_valid(self, data):
        pfCands = data["PFCandidate"]
        is_good = (
            (pfCands.pt > 1)
            & (abs(pfCands.eta) < 2.4)
            & (pfCands.hasTrackDetails)
        )
        pfCands_selected = pfCands[is_good]
        sort_idx = dask_awkward.argsort(pfCands_selected.pt, axis=-1, ascending=False)
        return pfCands_selected[sort_idx]
    
    def match_jet_to_pfcand(self, data, jet_name = None, pf_name = None, dR = 0.4):
        jets = data[jet_name]
        pfcands = data[pf_name]
        # here return_combinations=True is needed to return _pfcands_unzipped
        # which broadcasted in the way to duplicate pdcands list per every jet in event
        (dr, (_, _pfcands_unzipped)) = jets.metric_table(pfcands, metric=coffea.nanoevents.methods.vector.LorentzVector.delta_r, return_combinations=True)
        pfcands_matched = _pfcands_unzipped[(dr < dR)]
        return pfcands_matched
    
    def get_matched_pfCands(self, data, match_object, pf_name, dR = 0.4):
        pfCands = self.match_jet_to_pfcand(data, jet_name = match_object, pf_name = pf_name, dR = dR)
        pfCands_lead = awkward.firsts(pfCands, axis=-1)
        #pfCands_lead["dxysig"] = pfCands_lead.dxy / pfCands_lead.dxyError
        #pfCands_lead["Lrel"] = np.sqrt(pfCands_lead.dxy**2 + pfCands_lead.dz**2)
        #pfCands_lead["dxysig_weight"] = dask_awkward.mean(pfCands.dxy / pfCands.dxyError, weight=pfCands.pt, axis=-1)
        ## Sorting pfCands by dxy within each jet
        #sort_indices = dask_awkward.argsort(np.abs(pfCands.dxy), axis=-1, ascending=False)
        #sorted_pfCands = pfCands[sort_indices]
        #pfCands_leaddxy = awkward.firsts(sorted_pfCands, axis=-1)
        #pfCands_lead["maxdxysig"] = pfCands_leaddxy.dxy / pfCands_leaddxy.dxyError
        #pfCands_lead["maxdxy"] = pfCands_leaddxy.dxy
        return pfCands_lead
    
    def set_jet_dxy(self, jet_name, data):
        pf_name = "myPFCandidate"
        data[pf_name] = self.pfcand_valid(data = data)
        data["Jet_lead_pfcand"] = self.get_matched_pfCands(data = data, match_object = jet_name, pf_name = pf_name, dR=0.4)
        jets = data[jet_name]
        # Mask jets with dxy nan (no selected pfcands matching)
        bad_jets = dask_awkward.is_none(data["Jet_lead_pfcand"].dxy, axis=-1)
        jets = dask_awkward.mask(jets, ~bad_jets) # mask bad jets to keep coorect shape
        jets["dz"] = data["Jet_lead_pfcand"].dz
        jets["dxy"] = data["Jet_lead_pfcand"].dxy
        #jets["dxy_weight"] = np.abs(data["Jet_lead_pfcand"].dxy_weight)
        #jets["dxysig"] = np.abs(data["Jet_lead_pfcand"].dxysig)
        #jets["dxysig_weight"] = np.abs(data["Jet_lead_pfcand"].dxysig_weight)
        #jets["ip3d"] = np.sqrt(data["Jet_lead_pfcand"].dxy**2 + data["Jet_lead_pfcand"].dz**2)
        #jets["maxdxy"] = np.abs(data["Jet_lead_pfcand"].maxdxy)
        #jets["maxdxysig"] = np.abs(data["Jet_lead_pfcand"].maxdxysig)
        jets = jets[~bad_jets] # remove bad jets
        #jets = jets[jets.dxy >= self.config["jet_dxy_min"]]
        return jets
    
    
    # we will receive a NanoEvents instead of a coffea DataFrame
    def process(self, events) :
        
        output = self.accumulator
        
        sel_idx = (dask_awkward.num(events.GenVisTau, axis = 1) >= 2)
        events = events[sel_idx]
        
        jets = events.Jet[
            (events.Jet.pt > 30)
            & (abs(events.Jet.eta) < 2.4)
            & (events.Jet.jetId >= 4)
        ]
        
        events["Jet"] = jets
        events["Jet"] = self.set_jet_dxy(jet_name = "Jet", data = events)
        
        sel_idx = (dask_awkward.num(events.Jet, axis = 1) >= 2)
        events = events[sel_idx]
        
        events["Jet1"] = events.Jet[:, 0]
        events["Jet2"] = events.Jet[:, 1]
        
        # Skip processing as it is an EmptyArray
        if not dask_awkward.num(events, axis = 0) :
            
            return output
        
        print(self.dataset_args)
        dataset = events.metadata["dataset"]
        args = self.dataset_args.get(dataset, {})
        
        ctau0 = args.get("ctau0", None)
        ctau0_target = args.get("ctau0_target", None)
        
        print(args)
        print(ctau0, ctau0_target)
        
        do_lifetime_reweight = ctau0 and ctau0_target
        
        GenTau = events.GenPart[
            (abs(events.GenPart.pdgId) == 15)
            & events.GenPart.hasFlags(["isHardProcess"])
            & events.GenPart.hasFlags(["isFirstCopy"])
            #& ((abs(events.GenPart.distinctParent.pdgId) == 1000015) | (abs(events.GenPart.distinctParent.pdgId) == 2000015))
        ]
        events["GenTau"] = GenTau
        
        GenTau_sorted = events.GenTau[dask_awkward.argsort(GenTau.pt, axis = 1, ascending = False)]
        events["GenTau1"] = GenTau_sorted[:, 0]
        events["GenTau2"] = GenTau_sorted[:, 1]
        
        GenStau = GenTau.distinctParent
        events["GenStau"] = GenStau
        
        print("GenStau.pdgId:", GenStau.pdgId)
        print("num(events.GenStau):", dask_awkward.num(events.GenStau, axis = 1))
        
        events["GenStau", "gamma"] = events.GenStau.energy / events.GenStau.mass
        events["GenStau", "beta"] = (1 - 1/events.GenStau.gamma**2)**0.5
        events["GenStau", "pathL"] = ((GenStau.vertexX - events.GenTau.vertexX)**2 + (GenStau.vertexY - events.GenTau.vertexY)**2 + (GenStau.vertexZ - events.GenTau.vertexZ)**2)**0.5
        events["GenStau", "cpathL"] = events.GenStau.pathL / (events.GenStau.gamma * events.GenStau.beta)
        
        GenStau_sorted = events.GenStau[dask_awkward.argsort(GenStau.pt, axis = 1, ascending = False)]
        events["GenStau1"] = GenStau_sorted[:, 0]
        events["GenStau2"] = GenStau_sorted[:, 1]
        
        args_weight_perStau = {}
        args_weight_perEvent = {}
        
        if (do_lifetime_reweight) :
            
            #events["GenStau", "lifetime_weight"] = ctau0_target/ctau0 * numpy.exp(-events.GenStau.cpathL/ctau0_target) / numpy.exp(-events.GenStau.cpathL/ctau0)
            events["GenStau", "lifetime_weight"] = numpy.exp(-events.GenStau.cpathL/ctau0_target) / numpy.exp(-events.GenStau.cpathL/ctau0)
            events["lifetime_weight"] = dask_awkward.prod(events.GenStau.lifetime_weight, axis = 1)
            
            args_weight_perStau["weight"] = dask_awkward.flatten(events.GenStau.lifetime_weight)
            args_weight_perEvent["weight"] = events.lifetime_weight
        #
        #output["GenStau_pt"].fill(dataset = events.metadata["dataset"], GenStau_pt = dask_awkward.flatten(events.GenStau.pt), **args_weight_perStau)
        #output["GenStau_pathL"].fill(dataset = events.metadata["dataset"], GenStau_pathL = dask_awkward.flatten(events.GenStau.pathL), **args_weight_perStau)
        #output["GenStau_cpathL"].fill(dataset = events.metadata["dataset"], GenStau_cpathL = dask_awkward.flatten(events.GenStau.cpathL), **args_weight_perStau)
        
        for obj, d_qty in self.d_hist_scheme.items() :
            
            for qty, ax_key in d_qty.items() :
                
                d_ax_args = {}
                d_ax_args.update(self.d_hist_axis[ax_key])
                
                ax_name = f"{obj}_{qty}"
                
                d_hist_args = {"dataset": events.metadata["dataset"]}
                d_hist_args[ax_name] = events[obj][qty]
                d_hist_args.update(args_weight_perEvent)
                output[ax_name].fill(**d_hist_args)
        
        #output["GenStau"].fill(
        #    dataset = events.metadata["dataset"],
        #    GenStau_pt = dask_awkward.flatten(events.GenStau.pt),
        #    GenStau_pathL = dask_awkward.flatten(events.GenStau.pathL),
        #    GenStau_cpathL = dask_awkward.flatten(events.GenStau.cpathL),
        #    **args_weight_perStau
        #)
        #
        #output["GenStau1"].fill(
        #    dataset = events.metadata["dataset"],
        #    GenStau1_pt = events.GenStau1.pt,
        #    GenStau1_pathL = events.GenStau1.pathL,
        #    GenStau1_cpathL = events.GenStau1.cpathL,
        #    **args_weight_perEvent
        #)
        #
        #output["GenStau2"].fill(
        #    dataset = events.metadata["dataset"],
        #    GenStau2_pt = events.GenStau2.pt,
        #    GenStau2_pathL = events.GenStau2.pathL,
        #    GenStau2_cpathL = events.GenStau2.cpathL,
        #    **args_weight_perEvent
        #)
        #
        #output["GenTau"].fill(
        #    dataset = events.metadata["dataset"],
        #    GenTau_pt = dask_awkward.flatten(events.GenTau.pt),
        #    GenTau_vertexR = dask_awkward.flatten(events.GenTau.vertexR),
        #    **args_weight_perStau
        #)
        #
        #output["GenTau1"].fill(
        #    dataset = events.metadata["dataset"],
        #    GenTau1_pt = events.GenTau1.pt,
        #    GenTau1_vertexR = events.GenTau1.vertexR,
        #    **args_weight_perEvent
        #)
        #
        #output["GenTau2"].fill(
        #    dataset = events.metadata["dataset"],
        #    GenTau2_pt = events.GenTau2.pt,
        #    GenTau2_vertexR = events.GenTau2.vertexR,
        #    **args_weight_perEvent
        #)
        #
        #output["Jet1"].fill(
        #    dataset = events.metadata["dataset"],
        #    Jet1_pt = events.Jet1.pt,
        #    Jet1_disTauTag_score1 = events.Jet1.disTauTag_score1,
        #    **args_weight_perEvent
        #)
        #
        #output["Jet2"].fill(
        #    dataset = events.metadata["dataset"],
        #    Jet2_pt = events.Jet2.pt,
        #    Jet2_disTauTag_score1 = events.Jet2.disTauTag_score1,
        #    **args_weight_perEvent
        #)
        #
        #output["MET"].fill(
        #    dataset = events.metadata["dataset"],
        #    MET_pt = events.MET.pt,
        #    MET_px = events.MET.px,
        #    MET_py = events.MET.py,
        #    **args_weight_perEvent
        #)
        
        return output
    
    
    def postprocess(self, accumulator):
        
        pass
        #return accumulator





def main() :
    
    base_storage_dir = "/pnfs/desy.de/cms/tier2/store/user/sobhatta/LongLivedStaus/NanoAOD/Stau_2018UL"
    
    datasets = {
        #"MStau-250_ctau-1mm"  : f"{base_storage_dir}/SMS-TStauStau_MStau-250_ctau-1mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_250_ctau_1mm_mLSP_1/231218_144732/*/*.root",
        #"MStau-250_ctau-5mm"  : f"{base_storage_dir}/SMS-TStauStau_MStau-250_ctau-5mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_250_ctau_5mm_mLSP_1/231218_144739/*/*.root",
        #"MStau-250_ctau-10mm" : f"{base_storage_dir}/SMS-TStauStau_MStau-250_ctau-10mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_250_ctau_10mm_mLSP_1/231218_144745/*/*.root",
        "MStau-250_ctau-50mm" : f"{base_storage_dir}/SMS-TStauStau_MStau-250_ctau-50mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_250_ctau_50mm_mLSP_1/231218_144751/*/*.root",
        "MStau-250_ctau-100mm": f"{base_storage_dir}/SMS-TStauStau_MStau-250_ctau-100mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_250_ctau_100mm_mLSP_1/231218_144758/*/*.root",
        
        #"MStau-250_ctau-100mm": "nanoaod_with-disTauTagScore_1.root",
        
        #"MStau-250_ctau-1mm-from-10mm"      : f"{base_storage_dir}/SMS-TStauStau_MStau-250_ctau-10mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_250_ctau_10mm_mLSP_1/231218_144745/*/*.root",
        #"MStau-250_ctau-5mm-from-10mm"      : f"{base_storage_dir}/SMS-TStauStau_MStau-250_ctau-10mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_250_ctau_10mm_mLSP_1/231218_144745/*/*.root",
        #
        #"MStau-250_ctau-1mm-from-100mm"     : f"{base_storage_dir}/SMS-TStauStau_MStau-250_ctau-100mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_250_ctau_100mm_mLSP_1/231218_144758/*/*.root",
        #"MStau-250_ctau-5mm-from-100mm"     : f"{base_storage_dir}/SMS-TStauStau_MStau-250_ctau-100mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_250_ctau_100mm_mLSP_1/231218_144758/*/*.root",
        "MStau-250_ctau-50mm-from-100mm"    : f"{base_storage_dir}/SMS-TStauStau_MStau-250_ctau-100mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_250_ctau_100mm_mLSP_1/231218_144758/*/*.root",
    }
    
    dataset_args = {
        "MStau-250_ctau-1mm-from-10mm"      : {"ctau0_target": 0.1, "ctau0": 1.0},
        "MStau-250_ctau-5mm-from-10mm"      : {"ctau0_target": 0.5, "ctau0": 1.0},
        
        "MStau-250_ctau-1mm-from-100mm"     : {"ctau0_target": 0.1, "ctau0": 10.0},
        "MStau-250_ctau-5mm-from-100mm"     : {"ctau0_target": 0.5, "ctau0": 10.0},
        "MStau-250_ctau-50mm-from-100mm"    : {"ctau0_target": 5.0, "ctau0": 10.0},
    }
    
    for key, val in datasets.items() :
        
        l_file_tmp = glob.glob(val)[0: 1]
        datasets[key] = [f"{_ele}:Events" for _ele in l_file_tmp]
        #datasets[key] = [f"file://{_ele}:Events" for _ele in l_file_tmp]
    
    #print(datasets)
    
    #output = coffea.processor.run_uproot_job(
    #    datasets,
    #    "Events",
    #    MyProcessor(
    #        datasets = list(datasets.keys()),
    #        dataset_args = dataset_args,
    #    ),
    #    #executor = coffea.processor.iterative_executor,
    #    executor = coffea.processor.futures_executor,
    #    executor_args = {
    #        "schema": NanoAODSchema,
    #        #"skipbadfiles": True,
    #        #"xrootdtimeout": 600, #sec
    #        "workers": 5,
    #    },
    #    #chunksize=1000
    #)
    
    dataset_runnable, dataset_updated = coffea.dataset_tools.preprocess(
        fileset = datasets,
        align_clusters = False,
        step_size = 100_000,
        #files_per_batch = 1,
        skip_bad_files = True,
        save_form = False,
    )
    
    to_compute = coffea.dataset_tools.apply_to_fileset(
        data_manipulation = MyProcessor(
            #datasets = list(datasets.keys()),
            dataset_args = dataset_args,
        ),
        #coffea.dataset_tools.max_chunks(dataset_runnable, 300),
        fileset = dataset_runnable,
        schemaclass = NanoAODSchema,
    )
    
    ProgressBar().register()
    (output,) = dask.compute(to_compute)#, num_workers = 5)
    
    print(output)
    
    with uproot.recreate("output_stau_lifetime-reweighting_study.root") as fout:
        
        for dataset_key, dataset_output in output.items() :
            
            for hist_key, histo in dataset_output.items() :
                
                for ax in histo.axes :
                    
                    if ax.name == "dataset" :
                        continue
                    
                    print(dataset_key, ax)
                    fout[f"{dataset_key}/{ax.name}"] = histo[{"dataset": dataset_key}].project(ax.name)
    
    #with uproot.recreate("output_stau_lifetime-reweighting_study.root") as fout:
    #    
    #    for dset_key, dset_group in output.items():
    #        
    #        for hist_key, histo in dset_group.items() :
    #        #histo = output[ ]
    #        #print(key)
    #        #print(histo.axes)
    #        ##print(histo.axes["dataset"])
    #        #print(histo.__dict__)
    #        #print(histo.axes.__dict__)
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
        
    
    return 0


if __name__ == "__main__" :
    
    main()
