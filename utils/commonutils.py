#!/usr/bin/env python3

import json
import logging
import numpy
import re
import yaml
import ROOT


logging.basicConfig(format = "[%(levelname)s] [%(asctime)s] %(message)s", level = logging.INFO)
logger = logging.getLogger("mylogger")


# Factorize to a utils file later
def load_config(cfgfile) :
    
    #if (os.path.isfile(cfgfile)) :
    
    with open(cfgfile, "r") as fopen :
        
        content = fopen.read()
        
        if (cfgfile.endswith(".yml") or cfgfile.endswith(".yaml")) :
            
            logger.info(f"Loading yaml config: {cfgfile}")
            d_loadcfg = yaml.load(content, Loader = yaml.FullLoader)
        
        elif (cfgfile.endswith(".json")) :
            
            logger.info(f"Loading json config: {cfgfile}")
            d_loadcfg = json.loads(content)
        
        else :
            
            logger.error(f"Invalid config provided: {cfgfile}")
            exit(1)
    
    return d_loadcfg

def parse_string_regex(
    s,
    regexp,
    ) :
    
    rgx = re.compile(regexp)
    result = [m.groupdict() for m in rgx.finditer(s)][0]
    
    return result

def parse_stau_samplestring(
    s,
    regexp,
    ) :
    
    #rgx = re.compile("stau(\d+)_lsp(\d+)_ctau(\w+)")
    rgx = re.compile(regexp)
    
    #mstau, mlsp, ctau = rgx.findall(s)[0]
    #
    #result = {
    #    "mstau": float(mstau),
    #    "mlsp": float(mlsp),
    #    "ctau": ctau,
    #}
    
    result = [m.groupdict() for m in rgx.finditer(s)][0]
    
    # VERY hacky -- needs to be improved
    result["mstau"] = float(result["mstau"])
    result["mlsp"] = float(result["mlsp"])
    
    return result


def load_stau_xsec_file(xsecfile) :
    
    arr_xsec = numpy.loadtxt(xsecfile, delimiter = ",", converters = {1: eval})
    
    gr_xsec = ROOT.TGraph(len(arr_xsec[:, 0]), arr_xsec[:, 0], arr_xsec[:, 1])
    gr_xsec.SetName("gr_xsec")
    
    gr_xsec_unc_d = ROOT.TGraph(len(arr_xsec[:, 0]), arr_xsec[:, 0], arr_xsec[:, 2])
    gr_xsec_unc_d.SetName("gr_xsec_unc_d")
    
    gr_xsec_unc_u = ROOT.TGraph(len(arr_xsec[:, 0]), arr_xsec[:, 0], arr_xsec[:, 3])
    gr_xsec_unc_u.SetName("gr_xsec_unc_u")



def get_stau_xsec_dict(l_samplestr, xsecfile, regexp) :
    
    #arr_data = numpy.loadtxt(xsecfile, delimiter = ",", converters = {1: eval})
    arr_data = numpy.loadtxt(xsecfile, delimiter = ",")
    
    nelements = arr_data.shape[0]
    arr_mass = numpy.array(arr_data[:, 0], dtype = numpy.float64)
    arr_xsec = numpy.array(arr_data[:, 1], dtype = numpy.float64)
    
    #gr_xsec = ROOT.TGraph(nelements, arr_mass, arr_xsec)
    #gr_xsec.SetName("gr_xsec")
    
    #xsec = None
    #
    ## 1st column is the stau mass
    #rowidx = numpy.where(arr_xsec[:, 0] == mstau)[0]
    #
    ## Empty, i.e. mstau not present
    #if (not len(rowidx)) :
    #    
    #    logger.error(f"mstau {mstau} not found in {xsecfile}")
    #    exit(1)
    #
    #rowidx = rowidx[0]
    #
    ## 2nd colum is the xsec
    #xsec = arr_xsec[rowidx, 1]
    
    #xsec = gr_xsec.Eval(mstau)
    
    #d_xsec = {}
    #
    #for samplestr in l_samplestr :
    #    
    #    d_stau_param = parse_stau_samplestring(s = samplestr, regexp = regexp)
    #    mstau = d_stau_param["mstau"]
    #    
    #    fn = ROOT.TF1("xsec_fn", "[0] + ([1]*x**[2])", min(arr_mass), max(arr_mass))
    #    fit_result = gr_xsec.Fit(fn, "S")
    #    fitted_fn = gr_xsec.GetFunction("xsec_fn")
    #    
    #    #xsec = gr_xsec.Eval(mstau)
    #    xsec = fitted_fn.Eval(mstau)
    #    d_xsec[samplestr] = xsec
    #
    #return d_xsec
    
    d_xsec = {}
    
    for samplestr in l_samplestr :
        
        d_stau_param = parse_stau_samplestring(s = samplestr, regexp = regexp)
        mstau = d_stau_param["mstau"]
        
        # 1st column is the stau mass
        rowidx = numpy.where(arr_data[:, 0] == mstau)[0]
        
        # Empty, i.e. mstau not present
        if (not len(rowidx)) :
            
            logger.error(f"mstau {mstau} not found in {xsecfile}")
            exit(1)
        
        # 2nd colum is the xsec
        xsec = arr_data[rowidx, 1]
        
        d_xsec[samplestr] = xsec
    
    return d_xsec


#def get_stau_xsec_dict(l_samplestr, xsecfile, regexp) :
#    
#    d_xsec = {}
#    
#    l_mstaus = []
#    
#    for samplestr in l_samplestr :
#        
#        d_stau_param = parse_stau_samplestring(s = samplestr, regexp = regexp)
#        mstau = d_stau_param["mstau"]
#        l_mstaus.append(mstau)
#        
#        d_xsec[samplestr] = get_stau_xsec(samplestr = samplestr, xsecfile = xsecfile, regexp = regexp)
    
    return d_xsec


def get_hist(
    histfile,
    histname,
    samples,
    scales = [],
    rebin = None,
    underflow = True,
    overflow = True,
    set_min = None,
    set_max = None,
) :
    
    if (len(scales)) :
        
        assert(len(scales) == len(samples))
    
    hist_result = None
    
    for isample, sample in enumerate(samples) :
        
        histname_sample = f"{sample}/{histname}"
        logger.info(f"Getting [histogram {histname_sample}] from [file {histfile.GetName()}]")
        hist_sample = histfile.Get(histname_sample).Clone()
        hist_sample.SetDirectory(0)
        #hist_sample.Print("range")
        
        if (len(scales)) :
            
            hist_sample.Scale(scales[isample])
        
        if (rebin is not None) :
            
            rebin = numpy.array(rebin, dtype = numpy.float64)
            hist_sample = hist_sample.Rebin(len(rebin)-1, "", rebin)
        
        if (hist_result is None) :
            
            hist_result = hist_sample.Clone()
            hist_result.SetDirectory(0)
        
        else :
            
            hist_result.Add(hist_sample)
    
    
    # Add under/overflow if rebinned
    if (rebin is not None) :
        
        nbins = hist_result.GetNbinsX()
        
        if (underflow) :
            
            hist_result.AddBinContent(1, hist_result.GetBinContent(0))
            hist_result.SetBinContent(0, 0.0)
            hist_result.SetBinError(0, 0.0)
        
        if (overflow) :
            
            hist_result.AddBinContent(nbins, hist_result.GetBinContent(nbins+1))
            hist_result.SetBinContent(nbins+1, 0.0)
            hist_result.SetBinError(nbins+1, 0.0)
    
    if ((set_min is not None) or (set_max is not None)) :
        
        for ibin in range(hist_result.GetNbinsX()) :
            
            val = hist_result.GetBinContent(ibin+1)
            
            if ((set_min is not None) and (val < set_min)) :
                
                hist_result.SetBinContent(ibin+1, set_min)
            
            elif ((set_max is not None) and (val > set_max)) :
                
                hist_result.SetBinContent(ibin+1, set_max)
    
    return hist_result