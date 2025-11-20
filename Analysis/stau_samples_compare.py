#!/usr/bin/env python3

import awkward
import coffea
import coffea.processor
import dataclasses
import glob
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
            "pt"        : {"bins": 100, "start": 0, "stop": 1000},
            "pxy"       : {"bins": 200, "start": -1000, "stop": 1000},
            "pathL"     : {"bins": 1000, "start": 0, "stop": 100},
            "score"     : {"bins": 100, "start": 0, "stop": 1},
            "dxy"       : {"bins": 1000, "start": 0, "stop": 100},
            
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
        
        output = self.accumulator
        
        sel_idx = (awkward.num(events.GenVisTau, axis = 1) >= 2)
        events = events[sel_idx]
        
        jets = events.Jet[
            (events.Jet.pt > 30)
            & (abs(events.Jet.eta) < 2.4)
            & (events.Jet.jetId >= 4)
        ]
        
        events["Jet"] = jets
        events["Jet"] = self.set_jet_dxy(jet_name = "Jet", data = events)
        
        sel_idx = (awkward.num(events.Jet, axis = 1) >= 2)
        events = events[sel_idx]
        
        events["Jet1"] = events.Jet[:, 0]
        events["Jet2"] = events.Jet[:, 1]
        
        # Skip processing as it is an EmptyArray
        if not len(events) :
            
            return output
        
        #print(self.dataset_args)
        dataset = events.metadata["dataset"]
        args = self.dataset_args.get(dataset, {})
        
        ctau0 = args.get("ctau0", None)
        ctau0_target = args.get("ctau0_target", None)
        
        #print(args)
        #print(ctau0, ctau0_target)
        
        do_lifetime_reweight = ctau0 and ctau0_target
        
        GenTau = events.GenPart[
            (abs(events.GenPart.pdgId) == 15)
            & events.GenPart.hasFlags(["isHardProcess"])
            & events.GenPart.hasFlags(["isFirstCopy"])
            #& ((abs(events.GenPart.distinctParent.pdgId) == 1000015) | (abs(events.GenPart.distinctParent.pdgId) == 2000015))
        ]
        events["GenTau"] = GenTau
        
        GenTau_sorted = events.GenTau[awkward.argsort(GenTau.pt, axis = 1, ascending = False)]
        events["GenTau1"] = GenTau_sorted[:, 0]
        events["GenTau2"] = GenTau_sorted[:, 1]
        
        GenStau = GenTau.distinctParent
        events["GenStau"] = GenStau
        
        #print("GenStau.pdgId:", GenStau.pdgId)
        #print("num(events.GenStau):", awkward.num(events.GenStau, axis = 1))
        
        events["GenStau", "gamma"] = events.GenStau.energy / events.GenStau.mass
        events["GenStau", "beta"] = (1 - 1/events.GenStau.gamma**2)**0.5
        events["GenStau", "pathL"] = ((GenStau.vertexX - events.GenTau.vertexX)**2 + (GenStau.vertexY - events.GenTau.vertexY)**2 + (GenStau.vertexZ - events.GenTau.vertexZ)**2)**0.5
        events["GenStau", "cpathL"] = events.GenStau.pathL / (events.GenStau.gamma * events.GenStau.beta)
        
        GenStau_sorted = events.GenStau[awkward.argsort(GenStau.pt, axis = 1, ascending = False)]
        events["GenStau1"] = GenStau_sorted[:, 0]
        events["GenStau2"] = GenStau_sorted[:, 1]
        
        args_weight_perStau = {}
        args_weight_perEvent = {}
        
        if (do_lifetime_reweight) :
            
            #events["GenStau", "lifetime_weight"] = ctau0_target/ctau0 * numpy.exp(-events.GenStau.cpathL/ctau0_target) / numpy.exp(-events.GenStau.cpathL/ctau0)
            events["GenStau", "lifetime_weight"] = numpy.exp(-events.GenStau.cpathL/ctau0_target) / numpy.exp(-events.GenStau.cpathL/ctau0)
            events["lifetime_weight"] = awkward.prod(events.GenStau.lifetime_weight, axis = 1)
            
            args_weight_perStau["weight"] = awkward.flatten(events.GenStau.lifetime_weight)
            args_weight_perEvent["weight"] = events.lifetime_weight
        
        for obj, d_qty in self.d_hist_scheme.items() :
            
            for qty, ax_key in d_qty.items() :
                
                d_ax_args = {}
                d_ax_args.update(self.d_hist_axis[ax_key])
                
                ax_name = f"{obj}_{qty}"
                
                d_hist_args = {"dataset": events.metadata["dataset"]}
                d_hist_args[ax_name] = events[obj][qty]
                d_hist_args.update(args_weight_perEvent)
                output[ax_name].fill(**d_hist_args)
        
        return output
    
    
    def postprocess(self, accumulator):
        
        #pass
        return accumulator



def main() :
    
    d_fnamelist = {}
    d_fnamelist["stau250_lsp1_ctau100mm_ul18_private"] = glob.glob("/pnfs/desy.de/cms/tier2/store/user/sobhatta/LongLivedStaus/NanoAOD/SUS-RunIISummer20UL18GEN-stau250_lsp1_ctau100mm_v6/crab_stau250_lsp1_ctau100mm/230311_014339/*/nanoaod_with-disTauTagScore_*.root")
    d_fnamelist["stau250_lsp1_ctau100mm_ul18_central"] = glob.glob("/pnfs/desy.de/cms/tier2/store/user/sobhatta/LongLivedStaus/NanoAOD/Stau_2018UL/SMS-TStauStau_MStau-250_ctau-100mm_mLSP-1_TuneCP5_13TeV-madgraphMLM-pythia8/crab_SMS_TStauStau_MStau_250_ctau_100mm_mLSP_1/231218_144758/*/nanoaod_with-disTauTagScore_*.root")
    
    datasets = sortedcontainers.SortedDict({
        "stau250_lsp1_ctau100mm_ul18_private": d_fnamelist["stau250_lsp1_ctau100mm_ul18_private"][0: 200],
        "stau250_lsp1_ctau100mm_ul18_central": d_fnamelist["stau250_lsp1_ctau100mm_ul18_central"],
    })
    
    dataset_args = {
        "stau250_lsp1_ctau100mm_ul18_private": {},
        "stau250_lsp1_ctau100mm_ul18_central": {},
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
            "workers": 15,
        },
    )
    
    print(output)
    
    with uproot.recreate("output/stau_samples_compare/output_stau_samples_compare.root") as fout:
        
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
