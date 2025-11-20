#!/usr/bin/env python3


def main() :
    
    l_masses = [
        90,
        100,
        125,
        150,
        175,
        200,
        225,
        250,
        275,
        300,
        350,
        400,
        450,
        500,
        550,
        600,
    ]
    
    # (target, source)
    l_reweight_pairs = [
        (2, 5),
        (3, 5),
        (4, 5),
        
        (6.5, 10),
        (8.5, 10),
        
        (20, 50),
        (30, 50),
        (40, 50),
        
        (60, 100),
        (80, 100),
        
        (150, 1000),
        (200, 1000),
        (300, 1000),
        (400, 1000),
        (500, 1000),
        (700, 1000),
    ]
    
    l_samples = []
    
    for mass in l_masses :
        
        for ctau0_target, ctau0_source in l_reweight_pairs :
            
            ctau0_target = str(ctau0_target).replace(".", "p")
            ctau0_source = str(ctau0_source).replace(".", "p")
            
            sample = f"SMS-TStauStau_MStau-{mass}_ctau-{ctau0_target}mm_mLSP-1_from_ctau-{ctau0_source}mm"
            
            l_samples.append(sample)
    
    #print("\n".join(l_samples))
    
    print("For \"mc_datasets\" dict:")
    print("\n".join([f"\"{_smp}\": []," for _smp in l_samples]))
    
    print("\n")
    print("For \"sig_lifetime_reweighting_datasets\" list:")
    print("\n".join([f"\"{_smp}\"," for _smp in l_samples]))
    
    print(f"{len(l_samples)} samples")
    
    return 0

if __name__ == "__main__" :
    
    main()