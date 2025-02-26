#!/bin/bash

# Exit if a command fails
#set -Eeu -o pipefail
set -e

TYPE=llstau_maximally-mixed

DIR=configs/limits/${TYPE}

CARDS=(
    "2016_preVFP/config_datacard_BRT2.yaml"
    "2016_postVFP/config_datacard_BRT2.yaml"
    "2017/config_datacard_BRT2.yaml"
    "2018/config_datacard_BRT2.yaml"
)

#cp -av configs/limits/llstau_maximally-mixed/* configs/limits/llstau_mass-degenerate/
#
#for f in $(find configs/limits/llstau_mass-degenerate/ | grep "yaml$"); do
#    echo "Converting file: $f"
#    sed -i "s/maximally-mixed/mass-degenerate/g" $f
#done

replace(){
    findstr=$1
    replstr=$2
    fname=$3
    
    occurences=$(grep $findstr $fname | wc -l)
    
    if [ "$occurences" -eq "0" ]; then
        echo "Error: string ${findstr} not found in ${fname}"
        exit 1
    fi
    
    sed -i "s#${findstr}#${replstr}#g" $fname
}

for card in ${CARDS[@]}; do
    infile=${DIR}/${card}
    outfile=$(sed s/BRT2/BRT1/g <<< $infile)
    #echo $infile $outfile
    if [[ "$infile" == "$outfile" ]]; then
        echo "Error: input and output files are same"
        echo "Input: ${infile}"
        echo "Output: ${outfile}"
        exit 1
    fi
    
    cp -v $infile $outfile
    
    replace "BRT2" "BRT1" $outfile
    replace "/bin2/" "/bin1/" $outfile
    replace "bin1to2" "bin0to1" $outfile
    replace "sig_procs.yaml" "sig_procs_for-BRT1.yaml" $outfile
done

echo "Success!"