#!/bin/bash -x

ERA="$1"

SUFFIX="$2"

if [ -n "$SUFFIX" ]; then
    SUFFIX="_${SUFFIX}"
fi

#./plot_deltaNLL_DisTauSF.py \
#--input ~/nfs_dust/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/tmp/test_DisTauSF_mass_nbins1-60-80/DisTauSF/channels_all/eras_all/ZMT_wp-p*_dxy-gt-0p*/higgsCombine.DisTauSF.MultiDimFit.mH120.root \
#--output results/test_DisTauSF_mass_nbins1-60-80/channels_all/deltaNLL.root \
#--era ${ERA} \
#--title "1 #mu#tau_{h} bin, with #mu#mu CR"
#
#./plot_deltaNLL_DisTauSF.py \
#--input ~/nfs_dust/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/tmp/test_DisTauSF_mass_nbins2-60-70-80/DisTauSF/channels_all/eras_all/ZMT_wp-p*_dxy-gt-0p*/higgsCombine.DisTauSF.MultiDimFit.mH120.root \
#--output results/test_DisTauSF_mass_nbins2-60-70-80/channels_all/deltaNLL.root \
#--era ${ERA} \
#--title "2 #mu#tau_{h} bins, with #mu#mu CR"

./plot_deltaNLL_DisTauSF.py \
--input /home/soham/mnt/desy_dust/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/results/DisTauSF/DisTauSF_mass_nbins1-60-80${SUFFIX}/DisTauSF/channels_all/eras_${ERA}/ZMT_wp-p*_dxy-gt-0p*/scan/higgsCombine.DisTauSF.MultiDimFit.mH120.root \
--output results/DisTauSF/DisTauSF_mass_nbins1-60-80${SUFFIX}/DisTauSF/channels_all/eras_${ERA}/scan/deltaNLL.root \
--era ${ERA} \
--title "1 #mu#tau_{h} bin, without #mu#mu CR"

#./plot_deltaNLL_DisTauSF.py \
#--input ~/nfs_dust/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/tmp/test_DisTauSF_mass_nbins2-60-70-80_no-dy-rateparam/DisTauSF/channels_all/eras_all/ZMT_wp-p*_dxy-gt-0p*/higgsCombine.DisTauSF.MultiDimFit.mH120.root \
#--output results/test_DisTauSF_mass_nbins2-60-70-80_no-dy-rateparam/channels_all/deltaNLL.root \
#--era ${ERA} \
#--title "2 #mu#tau_{h} bins, without #mu#mu CR"
