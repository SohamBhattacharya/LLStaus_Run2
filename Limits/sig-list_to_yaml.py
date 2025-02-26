#!/usr/bin/env python3

import numpy
from ruamel.yaml import YAML


def main() :
    
    #sigfile = "configs/limits/sig_list_test.txt"
    #outfile = "configs/limits/sig_procs_test.yaml"
    
    #sigfile = "configs/limits/sig_list.txt"
    #outfile = "configs/limits/sig_procs.yaml"
    
    sigfile = "configs/distributions/sig_list_for-distributions.txt"
    outfile = "configs/distributions/sig_procs_for-distributions.yaml"
    
    l_sig = numpy.loadtxt(sigfile, dtype = str)
    #print(l_sig)
    
    d_sig = {}
    
    for sig in l_sig :
        
        sig = str(sig)
        
        d_sig[sig] = {
            "ismc": True,
            "issusy": True,
            "xsnorm": True,
            "scaleby": "1",
            "samples": [sig],
        }
    
    d_procs = {}
    d_procs["procs"] = d_sig
    
    yaml = YAML()
    yaml.preserve_quotes = True
    yaml.width = 1024
    
    #print(d_sig)
    
    with open(outfile, "w") as fopen:
        
        yaml.dump(d_procs, fopen)


if __name__ == "__main__" :
    
    main()
