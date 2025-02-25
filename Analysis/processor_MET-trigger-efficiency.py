
# This file illustrates how to implement a processor, realizing the selection
# steps and outputting histograms and a cutflow with efficiencies.
# Here we create a very simplified version of the ttbar-to-dilep processor.
# One can run this processor using
# 'python3 -m pepper.runproc --debug example_processor.py example_config.json'
# Above command probably will need a little bit of time before all cuts are
# applied once. This is because a chunk of events are processed simultaneously.
# You change adjust the number of events in a chunk and thereby the memory
# usage by using the --chunksize parameter (the default value is 500000).

from nis import match
from unittest import result
import pepper
import awkward as ak
from functools import partial
import numpy as np

from coffea.nanoevents import NanoAODSchema
# np.set_printoptions(threshold=np.inf)


# All processors should inherit from pepper.ProcessorBasicPhysics
class Processor(pepper.ProcessorBasicPhysics):
    # We use the ConfigTTbarLL instead of its base Config, to use some of its
    # predefined extras
    config_class = pepper.ConfigTTbarLL
    # schema_class = NanoAODSchema
    # schema_class.mixins.update({
    #     "PFCandidate": "PtEtaPhiMCollection",
    #     "SV": "PtEtaPhiMCollection",
    # })

    def __init__(self, config, eventdir):
        # Initialize the class, maybe overwrite some config variables and
        # load additional files if needed
        # Can set and modify configuration here as well
        config["histogram_format"] = "root"
        # Need to call parent init to make histograms and such ready
        super().__init__(config, eventdir)

        # It is not recommended to put anything as member variable into a
        # a Processor because the Processor instance is sent as raw bytes
        # between nodes when running on HTCondor.

    def zero_handler(func):
        def _function(self, data, *args, **kwargs):
            if len(data) > 0: return func(self, data, *args, **kwargs)
            else: return ak.Array([])
        return _function

    def process_selection(self, selector, dsname, is_mc, filler):
        # Implement the selection steps: add cuts, define objects and/or
        # compute event weights

        # Add a cut only allowing events according to the golden JSON
        # The good_lumimask method is specified in pepper.ProcessorBasicPhysics
        # It also requires a lumimask to be specified in config
        if not is_mc:
            selector.add_cut("Lumi", partial(
                self.good_lumimask, is_mc, dsname))
        
        if is_mc :
            if (
                dsname.startswith("WJetsToLNu") or
                dsname.startswith("W1JetsToLNu") or
                dsname.startswith("W2JetsToLNu") or
                dsname.startswith("W3JetsToLNu") or
                dsname.startswith("W4JetsToLNu") or
                dsname.startswith("WNJetsToLNu")
            ):
                selector.add_cut("Wjet_fixGenWeight", partial(self.fix_wjets_genweight))
                selector.add_cut("Wjet_reweighting", partial(self.do_w_jet_reweighting))
        
        
        selector.add_cut("MET_filters", partial(self.met_filters, is_mc))
        
        if "2018" in self.config["year"]:
            selector.add_cut("HEM_veto", partial(self.HEM_veto, is_mc = is_mc))
        
        # Only allow events that pass triggers specified in config
        # This also takes into account a trigger order to avoid triggering
        # the same event if it's in two different data datasets.
        pos_triggers_tag, neg_triggers_tag = pepper.misc.get_trigger_paths_for(
            dsname,
            is_mc,
            self.config["tag_trigger_map"],
            self.config["dataset_trigger_order"]
        )
        
        selector.add_cut("passed_tag_triggers", partial(
            self.passing_trigger,
            pos_triggers_tag,
            neg_triggers_tag
        ))

        # Select valid objects
        selector.set_column("trig_muons", partial(self.trig_muons)) 
        selector.set_column("muons_valid", partial(self.muons, l_match_coll = ["trig_muons"], match_dR = 0.1))
        #selector.set_column("jets_valid", partial(self.jets, l_clean_coll = ["muons_valid"], clean_dR = 0.4))
        selector.set_column("jets_valid", partial(self.jets, l_clean_coll = ["trig_muons"], clean_dR = 0.4))

        # Add basics cuts (at least 2 good jets)

        # Pick up leading jet (by score)
        # selector.set_column("jet_1", partial(self.leading_jet, order=0)) #1-st jet
        # selector.set_column("jet_2", partial(self.leading_jet, order=1)) #2-nd jet

        selector.add_cut("has_1_muon", partial(self.has_n_objs, collname = "muons_valid", min_count = 1))
        selector.add_cut("has_2_jets", partial(self.has_n_objs, collname = "jets_valid", min_count = 2))
        #selector.add_cut("passed_probe_triggers", self.passed_probe_triggers, probe_triggers = self.config["probe_triggers"])
        
        pos_triggers_probe, neg_triggers_probe = pepper.misc.get_trigger_paths_for(
            dsname,
            is_mc,
            self.config["probe_trigger_map"],
            self.config["dataset_trigger_order"]
        )
        
        selector.add_cut("passed_probe_triggers", partial(
            self.passing_trigger,
            pos_triggers_probe,
            neg_triggers_probe
        ))
    
    
    @zero_handler
    def fix_wjets_genweight(self, data):
        weights = ak.nan_to_num(1.0/data["genWeight"], nan = 0.0, posinf = 0.0, neginf = 0.0)
        return weights
    
    
    @zero_handler
    def do_w_jet_reweighting(self, data):
        njet = data["LHE"]["Njets"]
        #print(njet, type(njet))
        #print(self.config["W_jet_reweight"], type(self.config["W_jet_reweight"]))
        weights = ak.Array(self.config["W_jet_reweight"])[njet]
        #weights = self.config["W_jet_reweight"][njet]
        #print(weights)
        return weights
    
    
    @zero_handler
    def HEM_veto(self, data, is_mc):
        weight = np.ones(len(data), dtype=np.float32)
        jets = data["Jet"]
        electron = data["Electron"]
        electron_in15or16_hem = ( (electron.pt > 20) & (electron.eta > -3.0) & (electron.eta < -1.3) & (electron.phi > -1.57) & (electron.phi < -0.87) )
        jet_in15or16_hem = ( (jets.pt > 20) & (jets.eta > -3.2) & (jets.eta < -1.3) & (jets.phi > -1.77) & (jets.phi < -0.67) )
        in_hem = (ak.any(electron_in15or16_hem, axis=-1) | ak.any(jet_in15or16_hem, axis=-1))
        if is_mc:
            weight[in_hem] = (1-0.66)
        else:
            weight[in_hem] = 0.0
        return weight
    
    
    @zero_handler
    def match_and_clean(self, data, coll, l_match_coll = [], match_dR = 0.4, l_clean_coll = [], clean_dR = 0.4) :
        
        for collname in l_match_coll :
            
            coll = self.match_nearest(data, coll1 = coll, coll2 = collname, dR = match_dR, not_matched = False)
        
        for collname in l_clean_coll :
            
            coll = self.match_nearest(data, coll1 = coll, coll2 = collname, dR = clean_dR, not_matched = True)
        
        return coll
    
    
    @zero_handler
    def jets(self, data, l_match_coll = [], match_dR = 0.4, l_clean_coll = [], clean_dR = 0.4):
        jets = data["Jet"]

        is_good = (
            (jets.pt > 30)
            & (abs(jets.eta) < 2.4)
            & (jets.jetId >= 4)
        )
        
        jets = jets[is_good]
        
        result = self.match_and_clean(
            data = data,
            coll = jets,
            l_match_coll = l_match_coll,
            match_dR = match_dR,
            l_clean_coll = l_clean_coll,
            clean_dR = clean_dR,
        )
        
        return result
    
    
    @zero_handler
    def muons(self, data, l_match_coll = [], match_dR = 0.1, l_clean_coll = [], clean_dR = 0.1):
        muons = data["Muon"]
        
        is_good = (
            (muons.pt > 30)
            & (abs(muons.eta) < 2.4)
            & (muons.tightId)
        )
        
        muons =  muons[is_good]
        
        result = self.match_and_clean(
            data = data,
            coll = muons,
            l_match_coll = l_match_coll,
            match_dR = match_dR,
            l_clean_coll = l_clean_coll,
            clean_dR = clean_dR,
        )
        
        return result
    
    
    @zero_handler
    def trig_muons(self, data):
        trig_muons = data["TrigObj"]
        
        # Bit 3 should be on for 1mu
        # https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv11/2022postEE/doc_WZ_TuneCP5_13p6TeV_pythia8_Run3Summer22EENanoAODv11-126X_mcRun3_2022_realistic_postEE_v1-v1.html
        is_good = (
            (abs(trig_muons.id) == 13)
            #& (trig_muons.filterBits >= 8)
            & ((trig_muons.filterBits >> 3) >= 1)
        )
        #print(trig_muons.filterBits)
        #print(trig_muons.filterBits >> 3)
        #print((trig_muons.filterBits >> 3) & 1)
        #print((trig_muons.filterBits >> 3) >= 1)
        #print(trig_muons.filterBits & 8)
        #print(is_good)
        
        return trig_muons[is_good]
    
    
    @zero_handler
    def has_n_objs(self, data, collname, min_count):
        return ak.num(data[collname]) >= min_count
    
    
    #@zero_handler    
    #def passed_probe_triggers(self, data, probe_triggers):
    #    return data["MET"].pt > 100 #GeV
    
    
    @zero_handler
    def match_nearest(self, data, coll1=None, coll2=None, dR = 0.4, not_matched = False):
        obj1 = coll1
        obj2 = coll2
        
        if isinstance(coll1, str):
            obj1 = data[coll1]
        
        if isinstance(coll2, str):
            obj2 = data[coll2]
        
        matches, dRlist = obj1.nearest(obj2, return_metric=True, threshold=dR)
        if not_matched:
            idx_matches_jets = ak.is_none(matches, axis=1)
        else:
            idx_matches_jets = ~ak.is_none(matches, axis=1)
        # results = obj1.mask[idx_matches_jets] # this function convert False to None e.g: [[obj, obj, None, None,..],..]
        results = obj1[idx_matches_jets] # this function drops all False e.g: [[obj, obj,..],..]
        #sort_idx = ak.argsort(results.pt, axis=-1, ascending=False)
        #return results[sort_idx]
        
        return results
