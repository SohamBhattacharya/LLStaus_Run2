#!/usr/bin/env python3

import awkward
import coffea
import coffea.processor
#import coffea.dataset_tools
#import dask
#import dask_awkward
import dataclasses
import glob
import hist
import numpy
import os
import particle
import socket
import sortedcontainers
import typing
import uproot

import ROOT
ROOT.gROOT.SetBatch(True)

from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from dask.diagnostics import ProgressBar

import utils.commonutils as cmut


HOSTNAME = socket.gethostname()

@dataclasses.dataclass
class MyProcessor(coffea.processor.ProcessorABC) :
    
    #datasets         : typing.List[str]
    dataset_args     : typing.Dict
    
    def __post_init__(self) :
        
        self.dataset_axis = hist.axis.StrCategory([], growth = True, name = "dataset", label = "dataset")
        
        self.d_hist_axis = {
            "pt"        : {"bins": 100, "start": 0, "stop": 1000},
            "score"     : {"bins": 100, "start": 0, "stop": 1},
            "dxy"       : {"bins": 1000, "start": 0, "stop": 100},
            "flavor"    : {"bins": 1010, "start": -10, "stop": 1000},
        }
        
        self.d_hist_scheme = {
            "Jet": {"pt": "pt", "disTauTag_score1": "score", "dxy": "dxy", "flavor": "flavor"},
            
            #"Jet1": {"pt": "pt", "disTauTag_score1": "score", "dxy": "dxy"},
            #"Jet2": {"pt": "pt", "disTauTag_score1": "score", "dxy": "dxy"},
            
            "MET": {"pt": "pt"},
        }
        
        self.d_hist2d_scheme = {
            "Jet": [
                {"dxy": "dxy", "flavor": "flavor"},
            ],
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
                
                self._accumulator[ax_name] = hist.Hist(self.dataset_axis, hist.axis.Regular(**d_ax_args), name = ax_name)
        
        for obj, l_yx_pair in self.d_hist2d_scheme.items() :
            
            for d_yx_pair in l_yx_pair :
                
                ax_key_y, ax_key_x = list(d_yx_pair.keys())
                ax_name_y, ax_name_x = list(d_yx_pair.values())
                
                name = f"{obj}_{ax_key_y}_vs_{ax_key_x}"
                ax_name_x = f"{obj}_{ax_name_x}"
                ax_name_y = f"{obj}_{ax_name_y}"
                
                ax_x = hist.axis.Regular(
                    name = ax_name_x,
                    label = ax_name_x,
                    **self.d_hist_axis[ax_key_x]
                )
                
                ax_y = hist.axis.Regular(
                    name = ax_name_y,
                    label = ax_name_y,
                    **self.d_hist_axis[ax_key_y]
                )
                
                self._accumulator[name] = hist.Hist(self.dataset_axis, ax_x, ax_y, name = name)
        
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
    
    def get_charge(self, layout, **kwargs) :
        #return particle.Particle.from_pdgid(pdgid).charge
        #print(layout, type(layout), type(layout).__name__, len(layout))
        
        def eval_charge(pdgid) :
            
            charge = 0
            try :
                charge = particle.Particle.from_pdgid(pdgid).charge
            except Exception as e:
                pass
            return charge

        if hasattr(layout, "is_numpy") and layout.is_numpy:
            #print(layout.data)
            charges = awkward.contents.NumpyArray(
                numpy.array([eval_charge(_x) for _x in layout.data])
            )
            #print(charges)
            return charges
        
    
    def set_genpart_charge(self, data, genpart_name = "GenPart"):
        
        genparts = data[genpart_name]
        #print("pdgids:", pdgids, type(pdgids))
        #charges = numpy.array([particle.Particle.from_pdgid(pdgid).charge for pdgid in pdgids])
        genparts["charge"] = awkward.transform(
            transformation = self.get_charge,
            array = genparts.pdgId
        )
        
        return genparts
    
    def set_jet_flavor(self, jet_name, data):
        
        """
        1. Get lead PF candidate of each jet
        2. Get closest GenParticle (first copy) to the lead PF candidate
        3. Check the parent of the GenParticle
        """
        
        jets = data[jet_name]
        
        # add hadronic tau flavour:
        #tau_vis = data.GenVisTau[ ((data.GenVisTau.pt > 30) & (abs(data.GenVisTau.eta) < 2.4) &
        #                          (data.GenVisTau.parent.hasFlags(["fromHardProcess"])))
        #                        ]
        #matches_tauhad, _ = jets.nearest(tau_vis, return_metric=True, threshold=0.4)
        #matches_mu, _  = jets.nearest(data["gen_mu"], return_metric=True, threshold=0.4)
        #matches_ele, _ = jets.nearest(data["gen_ele"], return_metric=True, threshold=0.4)
        
        gen_part = data.GenPart[
            #(data.GenPart.hasFlags(["isLastCopy"]))
            (numpy.abs(data.GenPart.pdgId) > 8)
            & (numpy.abs(data.GenPart.pdgId) != 21)
            & (numpy.abs(data.GenPart.pdgId) != 22)
            & (numpy.abs(data.GenPart.pdgId) != 23)
            & (numpy.abs(data.GenPart.pdgId) != 24)
            
            & (numpy.abs(data.GenPart.pdgId) > 18)
            
            & (data.GenPart.charge != 0)
        ]
        #gen_part = data.GenPart
        #print("gen_part:", gen_part)
        #print("jets['lead_pfcand']:", jets["lead_pfcand"])

        matches_gen_part, _ = jets["lead_pfcand"].nearest(gen_part, return_metric=True, threshold=0.4)
        #matches_gen_part, _ = jets.nearest(gen_part, return_metric=True, threshold=0.4)
        #print("matches_gen_part:", matches_gen_part)

        #updFlavour = awkward.where(~awkward.is_none(matches_tauhad, axis=1), 15, jets.partonFlavour)
        #updFlavour = awkward.where(~awkward.is_none(matches_mu, axis=1), 13, updFlavour)
        #updFlavour = awkward.where(~awkward.is_none(matches_ele, axis=1), 11, updFlavour)
        
        flavor = awkward.where(~awkward.is_none(matches_gen_part, axis=1), numpy.abs(matches_gen_part.pdgId), 0)
        #print("flavor:", flavor)
        
        jets["flavor"] = flavor
        
        return jets
        
        #return flavor
    
    # we will receive a NanoEvents instead of a coffea DataFrame
    def process(self, events) :
        
        dataset = events.metadata["dataset"]
        output = self.accumulator
        
        output["nEvents"].fill(dataset = dataset, nEvents = numpy.ones(len(events)))
        
        jets = events.Jet[
            (events.Jet.pt > 30)
            & (abs(events.Jet.eta) < 2.1)
            & (events.Jet.jetId >= 6)
            & (events.Jet.btagDeepFlavB < 0.2783) #"ul2018": WpTuple(0.0490, 0.2783, 0.7100),
            & (events.Jet.disTauTag_score1 > 0.99)
        ]
        
        events["Jet"] = jets
        events["Jet"] = self.set_jet_dxy(jet_name = "Jet", data = events)
        
        events["Jet"] = events.Jet[
            (events.Jet.dxy > 0.2)
        ]
        
        sel_idx = (awkward.num(events.Jet, axis = 1) >= 1)
        events = events[sel_idx]
        
        # Skip processing as it is an EmptyArray
        if not len(events) :
            
            return output
        
        events["GenPart"] = self.set_genpart_charge(data = events, genpart_name = "GenPart")
        #print(events.GenPart.pdgId)
        #print(events.GenPart.charge)
        
        events["Jet"] = self.set_jet_flavor("Jet", events)
        
        #print("len(events):", len(events))
        #print("events.Jet:", events.Jet)
        #print("awkward.num(events.Jet, axis = 1):", awkward.num(events.Jet, axis = 1))
        #print("X"*10)
        
        #events["Jet1"] = events.Jet[:, 0]
        #events["Jet2"] = events.Jet[:, 1]
        #
        #print("events.GenJet:", events.GenJet)
        #print("events.Jet1.genJetIdx:", events.Jet1.genJetIdx)
        #print("events.Jet2.genJetIdx:", events.Jet2.genJetIdx)
        #
        #events["Jet1_GenJet"] = events.GenJet[events.Jet1.genJetIdx]
        #events["Jet2_GenJet"] = events.GenJet[events.Jet2.genJetIdx]
        
        #print(self.dataset_args)
        
        args = self.dataset_args.get(dataset, {})
        
        for obj, d_qty in self.d_hist_scheme.items() :
            
            for qty, ax_key in d_qty.items() :
                
                d_ax_args = {}
                d_ax_args.update(self.d_hist_axis[ax_key])
                
                ax_name = f"{obj}_{qty}"
                
                d_hist_args = {"dataset": events.metadata["dataset"]}
                d_hist_args[ax_name] = awkward.flatten(events[obj][qty], axis = None)
                output[ax_name].fill(**d_hist_args)
        
        for obj, l_yx_pair in self.d_hist2d_scheme.items() :
            
            for d_yx_pair in l_yx_pair :
                
                ax_key_y, ax_key_x = list(d_yx_pair.keys())
                ax_name_y, ax_name_x = list(d_yx_pair.values())
                
                qty_x = ax_key_x
                qty_y = ax_key_y
                
                name = f"{obj}_{ax_key_y}_vs_{ax_key_x}"
                ax_name_x = f"{obj}_{ax_name_x}"
                ax_name_y = f"{obj}_{ax_name_y}"
                
                d_hist_args = {"dataset": events.metadata["dataset"]}
                d_hist_args[ax_name_x] = awkward.flatten(events[obj][qty_x], axis = None)
                d_hist_args[ax_name_y] = awkward.flatten(events[obj][qty_y], axis = None)
                
                #print("d_hist_args:", d_hist_args)
                
                output[name].fill(**d_hist_args)
        
        return output
    
    
    def postprocess(self, accumulator):
        
        #pass
        return accumulator



def main() :
    
    d_fnamelist = {}
    
    l_paths = [
        #"/pnfs/desy.de/cms/tier2/store/user/sobhatta/LongLivedStaus/NanoAOD/QCD_Pt_15to30_TuneCP5_13TeV_pythia8/crab_QCD_Pt_15to30/230311_140730",
        
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
    
    if "naf" in HOSTNAME :
        
        for path in l_paths :
            
            sample_name = path.split("/")[-3]
            d_fnamelist[sample_name] = cmut.natural_sort(glob.glob(f"{path}/*/nanoaod_with-disTauTagScore_*.root"))#[0: 20]
        #d_fnamelist["QCD"] = cmut.natural_sort(glob.glob("/pnfs/desy.de/cms/tier2/store/user/sobhatta/LongLivedStaus/NanoAOD/QCD_Pt_50to80_TuneCP5_13TeV_pythia8/crab_QCD_Pt_50to80/230311_140741/*/nanoaod_with-disTauTagScore_*.root"))
    
    else :
        #d_fnamelist["QCD"] = cmut.natural_sort(glob.glob("nanoaod_SUS-RunIISummer20UL18GEN-stau250_lsp1_ctau100mm_v6_sobhatta-MiniAOD-c15273f0b6812ff053a850f456209388_USER_with-disTauTagScore.root"))
        d_fnamelist["QCD"] = cmut.natural_sort(glob.glob("data/NanoAOD/QCD_Pt_50to80_TuneCP5_13TeV_pythia8/crab_QCD_Pt_50to80/230311_140741/*/nanoaod_with-disTauTagScore_*.root"))
    
    datasets = sortedcontainers.SortedDict({
        #"QCD": d_fnamelist["QCD"][0: min(1000, len(d_fnamelist["QCD"]))],
        #"QCD1": d_fnamelist["QCD"][0: 5],

        _key: _val for _key, _val in d_fnamelist.items()
    })
    
    #datasets = sortedcontainers.SortedDict({
    #    "QCD": [f"file:{_f}:Events" for _f in d_fnamelist["QCD"][0: 1]],
    #})
    
    dataset_args = {
        _key: {} for _key in d_fnamelist.keys()
    }
    
    schema_class = NanoAODSchema
    schema_class.mixins.update({
        "PFCandidate": "PtEtaPhiMCollection",
        "SV": "PtEtaPhiMCollection",
        #"GenJet": "PtEtaPhiMCollection",
    })
    
    #if "naf" in HOSTNAME :
    #    
    #    output = coffea.processor.run_uproot_job(
    #        datasets,
    #        "Events",
    #        MyProcessor(
    #            dataset_args = dataset_args,
    #        ),
    #        executor = coffea.processor.futures_executor,
    #        executor_args = {
    #            "schema": NanoAODSchema,
    #            "skipbadfiles": True,
    #            "workers": 15,
    #        },
    #    )
    #
    #else :
    #    
    #    run = coffea.processor.Runner(
    #        executor = coffea.processor.FuturesExecutor(compression = None, workers = 15),
    #        schema = schema_class,
    #        chunksize = 10_000,
    #        #chunksize = 100,
    #        # maxchunks=10,  # total 676 chunks
    #        skipbadfiles = True
    #    )
    #    
    #    output = run(
    #        fileset = datasets,
    #        treename = "Events",
    #        processor_instance = MyProcessor(
    #            dataset_args = dataset_args,
    #        ),
    #    )
    
    run = coffea.processor.Runner(
        executor = coffea.processor.FuturesExecutor(compression = None, workers = 20),
        schema = schema_class,
        chunksize = 10_000,
        #chunksize = 100,
        # maxchunks=10,  # total 676 chunks
        skipbadfiles = True
    )
    
    output = run(
        fileset = datasets,
        treename = "Events",
        processor_instance = MyProcessor(
            dataset_args = dataset_args,
        ),
    )
    
    print("output:", output)
    #outfname = "output/fake-jets_study/fake-jets_study.root"
    #outfname = "output/fake-jets_study_w-charged-genparticles/fake-jets_study.root"
    #outfname = "output/fake-jets_study_w-charged-genparticles_w-bveto/fake-jets_study.root"
    outfname = "output/fake-jets_study_w-charged-genparticles_wo-genleptons_w-bveto/fake-jets_study.root"
    outdir = os.path.dirname(outfname)
    
    if outdir :
        os.system(f"mkdir -p {outdir}")
    
    with uproot.recreate(outfname) as fout:
            
        for dataset_key in datasets.keys() :
            
            for hist_key, histo in output.items() :
                
                print(dataset_key, histo.name, histo)
                
                if "_vs_" in histo.name :
                    
                    print(dataset_key, histo)
                    h_tmp = histo[{"dataset": dataset_key}]
                    fout[f"{dataset_key}/{histo.name}"] = h_tmp

                else :
                    for ax in histo.axes :
                        
                        if ax.name == "dataset" :
                            continue
                        
                        print(dataset_key, ax, histo)
                        h_tmp = histo[{"dataset": dataset_key}].project(ax.name)
                        fout[f"{dataset_key}/{ax.name}"] = h_tmp
    
    print(f"Output written to: {outfname}")
    
    return 0


if __name__ == "__main__" :
    
    main()
