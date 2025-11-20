#!/usr/bin/env python3

import awkward
import coffea
import coffea.processor
import dataclasses
import glob
from hepunits import d
import hist
import numpy
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
        
        self.dataset_axis = hist.axis.StrCategory([], growth = True, name = "dataset", label = "dataset")
        
        self.d_hist_axis = {
            "e": {"bins": 100, "start": 0, "stop": 1000},
            "pt": {"bins": 100, "start": 0, "stop": 1000},
            "eta": {"bins": 600, "start": -3.0, "stop": 3.0},
            "score": {"bins": 100, "start": 0, "stop": 1.0},
            "dxy": {"bins": 1000, "start": 0, "stop": 100},
        }
        
        # object: {quantity: axis_key}
        self.d_hist_scheme = {
            "Jet":          {"pt": "pt", "eta": "eta", "disTauTag_score1": "score", "dxy": "dxy"},
            "Jet_light":    {"pt": "pt", "eta": "eta", "disTauTag_score1": "score", "dxy": "dxy"},
            "Jet_b":        {"pt": "pt", "eta": "eta", "disTauTag_score1": "score", "dxy": "dxy"},
            "Jet_c":        {"pt": "pt", "eta": "eta", "disTauTag_score1": "score", "dxy": "dxy"},
            "Jet_bc":       {"pt": "pt", "eta": "eta", "disTauTag_score1": "score", "dxy": "dxy"},
            #"Jet1": {"pt": "pt"},
            #"Jet2": {"pt": "pt"},
            
            "MET":          {"pt": "pt"},
        }
        
        self._accumulator = {}
        
        self._accumulator["nEvents"] = hist.Hist(
            self.dataset_axis,
            hist.axis.Regular(name = "nEvents", label = "nEvents", bins = 1, start = 0, stop = 2),
            name = "nEvents"
        )
        
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
    
    
    def pfcand_valid(self, data):
        pfCands = data["PFCandidate"]
        is_good = (
            (pfCands.pt > 1)
            & (abs(pfCands.eta) < 2.4)
            & (pfCands.hasTrackDetails)
        )
        pfCands_selected = pfCands[is_good]
        sort_idx = awkward.argsort(pfCands_selected.pt, axis=-1, ascending=False)
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
        pfCands_lead["dxysig"] = pfCands_lead.dxy / pfCands_lead.dxyError
        #pfCands_lead["Lrel"] = numpy.sqrt(pfCands_lead.dxy**2 + pfCands_lead.dz**2)
        #pfCands_lead["dxysig_weight"] = awkward.mean(pfCands.dxy / pfCands.dxyError, weight=pfCands.pt, axis=-1)
        ## Sorting pfCands by dxy within each jet
        #sort_indices = awkward.argsort(numpy.abs(pfCands.dxy), axis=-1, ascending=False)
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
        jets["lead_pfcand"] = data["Jet_lead_pfcand"]
        # Mask jets with dxy nan (no selected pfcands matching)
        bad_jets = awkward.is_none(data["Jet_lead_pfcand"].dxy, axis=-1)
        jets = awkward.mask(jets, ~bad_jets) # mask bad jets to keep coorect shape
        jets["dz"] = numpy.abs(data["Jet_lead_pfcand"].dz)
        jets["dxy"] = numpy.abs(data["Jet_lead_pfcand"].dxy)
        #jets["dxy_weight"] = numpy.abs(data["Jet_lead_pfcand"].dxy_weight)
        jets["dxysig"] = numpy.abs(data["Jet_lead_pfcand"].dxysig)
        #jets["dxysig_weight"] = numpy.abs(data["Jet_lead_pfcand"].dxysig_weight)
        #jets["ip3d"] = numpy.sqrt(data["Jet_lead_pfcand"].dxy**2 + data["Jet_lead_pfcand"].dz**2)
        #jets["maxdxy"] = numpy.abs(data["Jet_lead_pfcand"].maxdxy)
        #jets["maxdxysig"] = numpy.abs(data["Jet_lead_pfcand"].maxdxysig)
        jets = jets[~bad_jets] # remove bad jets
        #jets = jets[jets.dxy >= self.config["jet_dxy_min"]]
        return jets
    
    # we will receive a NanoEvents instead of a coffea DataFrame
    def process(self, events) :
        
        dataset = events.metadata["dataset"]
        args = self.dataset_args.get(dataset, {})
        
        is_signal = args.get("is_signal", False)
        
        output = self.accumulator
        
        output["nEvents"].fill(dataset = dataset, nEvents = numpy.ones(len(events)))
        
        sel_idx = (awkward.num(events.GenVisTau, axis = 1) > 0)
        events = events[sel_idx]
        
        jets = events.Jet[
            (events.Jet.pt > 30)
            & (abs(events.Jet.eta) < 2.1)
            & (events.Jet.jetId >= 6)
            #& (events.Jet.dxy > 0.2)
        ]
        
        events["Jet"] = self.set_jet_dxy(jet_name = "Jet", data = events)
        jets = events.Jet
        #jets = events.Jet[(events.Jet.dxy > 0.2)]
        
        genObjs = events.GenVisTau if is_signal else events.GenJet
        
        (dr, (_, jets)) = genObjs.metric_table(jets, metric = coffea.nanoevents.methods.vector.LorentzVector.delta_r, return_combinations = True)
        jets = jets[(dr < 0.4)]
        jets = jets[~awkward.is_none(jets, axis = 1)]
        
        #print(jets)
        
        events["Jet"] = jets
        sel_idx = (awkward.num(events.Jet, axis = 1) > 0)
        events = events[sel_idx]
        
        events["Jet_light"] = jets[
            (jets.partonFlavour != 0)
            & (abs(jets.partonFlavour) != 4)
            & (abs(jets.partonFlavour) != 5)
        ]
        
        sel_c = (abs(jets.partonFlavour) == 4)
        sel_b = (abs(jets.partonFlavour) == 5)
        events["Jet_c"] = jets[sel_c]
        events["Jet_b"] = jets[sel_b]
        events["Jet_bc"] = jets[sel_b | sel_c]
        
        #events["Jet1"] = events.Jet[:, 0]
        #events["Jet2"] = events.Jet[:, 1]
        
        # Skip processing as it is an EmptyArray
        if not len(events) :
            
            return output
        
        for obj, d_qty in self.d_hist_scheme.items() :
            
            for qty, ax_key in d_qty.items() :
                
                d_ax_args = {}
                d_ax_args.update(self.d_hist_axis[ax_key])
                
                ax_name = f"{obj}_{qty}"
                
                d_hist_args = {"dataset": events.metadata["dataset"]}
                
                d_hist_args[ax_name] = awkward.flatten([events[obj][qty]], axis = None)
                
                if len(d_hist_args[ax_name]) :
                    
                    output[ax_name].fill(**d_hist_args)
        
        return output
    
    
    def postprocess(self, accumulator):
        
        return accumulator





def main() :
    
    
    base_storage_dir = "/pnfs/desy.de/cms/tier2/store/user/sobhatta/LongLivedStaus/NanoAOD"
    
    datasets = {
        
        #"MStau-100_ctau-1mm": f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-100_ctau-1mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_100_ctau_1mm_mLSP_1/231218_144035/*/*.root",
        #"MStau-100_ctau-5mm": f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-100_ctau-5mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_100_ctau_5mm_mLSP_1/231218_144041/*/*.root",
        #"MStau-100_ctau-10mm": f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-100_ctau-10mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_100_ctau_10mm_mLSP_1/231218_144048/*/*.root",
        #"MStau-100_ctau-50mm": f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-100_ctau-50mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_100_ctau_50mm_mLSP_1/231218_144055/*/*.root",
        #"MStau-100_ctau-100mm": f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-100_ctau-100mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_100_ctau_100mm_mLSP_1/231218_144102/*/*.root",
        #"MStau-100_ctau-1000mm": f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-100_ctau-1000mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_100_ctau_1000mm_mLSP_1/231218_144108/*/*.root",
        
        #"MStau-200_ctau-1000mm": f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-200_ctau-1000mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_200_ctau_1000mm_mLSP_1/231218_144548/*/*.root",
        "MStau-200_ctau-100mm": f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-200_ctau-100mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_200_ctau_100mm_mLSP_1/231218_144542/*/*.root",
        #"MStau-200_ctau-10mm": f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-200_ctau-10mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_200_ctau_10mm_mLSP_1/231218_144529/*/*.root",
        #"MStau-200_ctau-1mm": f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-200_ctau-1mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_200_ctau_1mm_mLSP_1/231218_144517/*/*.root",
        #"MStau-200_ctau-50mm": f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-200_ctau-50mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_200_ctau_50mm_mLSP_1/231218_144536/*/*.root",
        "MStau-200_ctau-5mm": f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-200_ctau-5mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_200_ctau_5mm_mLSP_1/231218_144523/*/*.root",
        
        #"MStau-250_ctau-1mm"  : f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-250_ctau-1mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_250_ctau_1mm_mLSP_1/231218_144732/*/*.root",
        #"MStau-250_ctau-5mm"  : f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-250_ctau-5mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_250_ctau_5mm_mLSP_1/231218_144739/*/*.root",
        #"MStau-250_ctau-10mm" : f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-250_ctau-10mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_250_ctau_10mm_mLSP_1/231218_144745/*/*.root",
        #"MStau-250_ctau-50mm" : f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-250_ctau-50mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_250_ctau_50mm_mLSP_1/231218_144751/*/*.root",
        #"MStau-250_ctau-100mm": f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-250_ctau-100mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_250_ctau_100mm_mLSP_1/231218_144758/*/*.root",
        #"MStau-250_ctau-1000mm": f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-250_ctau-1000mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_250_ctau_1000mm_mLSP_1/231218_144804/*/*.root",
        
        #"MStau-400_ctau-1mm": f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-400_ctau-1mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_400_ctau_1mm_mLSP_1/231218_145205/*/*.root",
        "MStau-400_ctau-5mm": f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-400_ctau-5mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_400_ctau_5mm_mLSP_1/231218_145213/*/*.root",
        #"MStau-400_ctau-10mm": f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-400_ctau-10mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_400_ctau_10mm_mLSP_1/231218_153817/*/*.root",
        #"MStau-400_ctau-50mm": f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-400_ctau-50mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_400_ctau_50mm_mLSP_1/231218_145227/*/*.root",
        "MStau-400_ctau-100mm": f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-400_ctau-100mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_400_ctau_100mm_mLSP_1/231218_153857/*/*.root",
        #"MStau-400_ctau-1000mm": f"{base_storage_dir}/Stau_2018UL/SMS-TStauStau_MStau-400_ctau-1000mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_400_ctau_1000mm_mLSP_1/231218_145240/*/*.root",
        
        "TTToHadronic" : f"{base_storage_dir}/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/crab_TTToHadronic/230311_084410/*/*.root",
    }
    
    l_paths = [
        "/pnfs/desy.de/cms/tier2/store/user/sobhatta/LongLivedStaus/NanoAOD/QCD_Pt_30to50_TuneCP5_13TeV_pythia8/crab_QCD_Pt_30to50/230311_140736",
        "/pnfs/desy.de/cms/tier2/store/user/sobhatta/LongLivedStaus/NanoAOD/QCD_Pt_50to80_TuneCP5_13TeV_pythia8/crab_QCD_Pt_50to80/230311_140741",
        "/pnfs/desy.de/cms/tier2/store/user/sobhatta/LongLivedStaus/NanoAOD/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/crab_QCD_Pt_80to120/230311_140747",
        "/pnfs/desy.de/cms/tier2/store/user/sobhatta/LongLivedStaus/NanoAOD/QCD_Pt_120to170_TuneCP5_13TeV_pythia8/crab_QCD_Pt_120to170/230311_140752",
        "/pnfs/desy.de/cms/tier2/store/user/sobhatta/LongLivedStaus/NanoAOD/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/crab_QCD_Pt_170to300/230311_140757",
        "/pnfs/desy.de/cms/tier2/store/user/sobhatta/LongLivedStaus/NanoAOD/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/crab_QCD_Pt_300to470/230311_140802",
        "/pnfs/desy.de/cms/tier2/store/user/sobhatta/LongLivedStaus/NanoAOD/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/crab_QCD_Pt_470to600/230311_140807",
        "/pnfs/desy.de/cms/tier2/store/user/sobhatta/LongLivedStaus/NanoAOD/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/crab_QCD_Pt_600to800/230311_140814",
        "/pnfs/desy.de/cms/tier2/store/user/sobhatta/LongLivedStaus/NanoAOD/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/crab_QCD_Pt_800to1000/230311_140821",
        "/pnfs/desy.de/cms/tier2/store/user/sobhatta/LongLivedStaus/NanoAOD/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/crab_QCD_Pt_1000to1400/230311_140827",
        "/pnfs/desy.de/cms/tier2/store/user/sobhatta/LongLivedStaus/NanoAOD/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/crab_QCD_Pt_1400to1800/230311_140832",
        "/pnfs/desy.de/cms/tier2/store/user/sobhatta/LongLivedStaus/NanoAOD/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/crab_QCD_Pt_1800to2400/230311_140838",
        "/pnfs/desy.de/cms/tier2/store/user/sobhatta/LongLivedStaus/NanoAOD/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/crab_QCD_Pt_2400to3200/230311_140844",
        "/pnfs/desy.de/cms/tier2/store/user/sobhatta/LongLivedStaus/NanoAOD/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/crab_QCD_Pt_3200toInf/230311_140850",
    ]
    
    for key, val in datasets.items() :
        
        l_file_tmp = glob.glob(val)
        l_file_tmp = l_file_tmp[0: min(200, len(l_file_tmp))]
        #datasets[key] = [f"{_ele}" for _ele in l_file_tmp]
        datasets[key] = l_file_tmp
    
    for path in l_paths :
        
        sample_name = path.split("/")[-3]
        l_file_tmp = glob.glob(f"{path}/*/nanoaod_with-disTauTagScore_*.root")
        l_file_tmp = l_file_tmp[0: min(200, len(l_file_tmp))]
        datasets[sample_name] = l_file_tmp
    
    dataset_args = {
        _key : {"is_signal": True if "Stau" in _key else False} for _key in datasets.keys()
    }
    
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
            "workers": 20,
        },
    )
    
    print(output)
    
    outdir = "output/stau_score_study"
    os.system(f"mkdir -p {outdir}")
    
    with uproot.recreate(f"{outdir}/output_stau_score_study.root") as fout:
        
        for dataset_key in datasets.keys() :
            
            for hist_key, histo in output.items() :
                
                for ax in histo.axes :
                    
                    if ax.name == "dataset" :
                        continue
                    
                    print(dataset_key, ax)
                    h_tmp = histo[{"dataset": dataset_key}].project(ax.name)
                    #h_tmp = h_tmp / h_tmp.sum(flow = True)
                    fout[f"{dataset_key}/{ax.name}"] = h_tmp
    
    return 0


if __name__ == "__main__" :
    
    main()
