#!/bin/bash

./plot_deltaNLL_DisTauSFs.py \
--input ~/nfs_dust/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/tmp/test_DisTauSF_mass_nbins1-60-80/DisTauSF/channels_all/eras_all/ZMT_wp-p*_dxy-gt-0p*/higgsCombine.DisTauSF.MultiDimFit.mH120.root \
--output results/test_DisTauSF_mass_nbins1-60-80/channels_all/deltaNLL.root \
--era 2018 \
--title "1 #mu#tau_{h} bin, with #mu#mu CR"

./plot_deltaNLL_DisTauSFs.py \
--input ~/nfs_dust/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/tmp/test_DisTauSF_mass_nbins2-60-70-80/DisTauSF/channels_all/eras_all/ZMT_wp-p*_dxy-gt-0p*/higgsCombine.DisTauSF.MultiDimFit.mH120.root \
--output results/test_DisTauSF_mass_nbins2-60-70-80/channels_all/deltaNLL.root \
--era 2018 \
--title "2 #mu#tau_{h} bins, with #mu#mu CR"

./plot_deltaNLL_DisTauSFs.py \
--input ~/nfs_dust/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/tmp/test_DisTauSF_mass_nbins1-60-80_no-dy-rateparam/DisTauSF/channels_all/eras_all/ZMT_wp-p*_dxy-gt-0p*/higgsCombine.DisTauSF.MultiDimFit.mH120.root \
--output results/test_DisTauSF_mass_nbins1-60-80_no-dy-rateparam/channels_all/deltaNLL.root \
--era 2018 \
--title "1 #mu#tau_{h} bin, without #mu#mu CR"

./plot_deltaNLL_DisTauSFs.py \
--input ~/nfs_dust/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/tmp/test_DisTauSF_mass_nbins2-60-70-80_no-dy-rateparam/DisTauSF/channels_all/eras_all/ZMT_wp-p*_dxy-gt-0p*/higgsCombine.DisTauSF.MultiDimFit.mH120.root \
--output results/test_DisTauSF_mass_nbins2-60-70-80_no-dy-rateparam/channels_all/deltaNLL.root \
--era 2018 \
--title "2 #mu#tau_{h} bins, without #mu#mu CR"
