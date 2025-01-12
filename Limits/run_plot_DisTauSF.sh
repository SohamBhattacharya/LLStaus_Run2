#!/bin/bash

./plot_DisTauSFs.py \
--jsons ~/nfs_dust/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/tmp/test_DisTauSF_mass_nbins1-60-80/DisTauSF/channels_all/eras_all/fit_result_ZMT_wp-p*.json \
--output results/test_DisTauSF_mass_nbins1-60-80/channels_all/DisTauSF.root \
--era 2018 \
--title "1 #mu#tau_{h} bin, with #mu#mu CR"


./plot_DisTauSFs.py \
--jsons ~/nfs_dust/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/tmp/test_DisTauSF_mass_nbins2-60-70-80/DisTauSF/channels_all/eras_all/fit_result_ZMT_wp-p*.json \
--output results/test_DisTauSF_mass_nbins2-60-70-80/channels_all/DisTauSF.root \
--era 2018 \
--title "2 #mu#tau_{h} bins, with #mu#mu CR"


./plot_DisTauSFs.py \
--jsons ~/nfs_dust/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/tmp/test_DisTauSF_mass_nbins1-60-80_no-dy-rateparam/DisTauSF/channels_all/eras_all/fit_result_ZMT_wp-p*.json \
--output results/test_DisTauSF_mass_nbins1-60-80_no-dy-rateparam/channels_all/DisTauSF.root \
--era 2018 \
--title "1 #mu#tau_{h} bin, without #mu#mu CR"


./plot_DisTauSFs.py \
--jsons ~/nfs_dust/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/tmp/test_DisTauSF_mass_nbins2-60-70-80_no-dy-rateparam/DisTauSF/channels_all/eras_all/fit_result_ZMT_wp-p*.json \
--output results/test_DisTauSF_mass_nbins2-60-70-80_no-dy-rateparam/channels_all/DisTauSF.root \
--era 2018 \
--title "2 #mu#tau_{h} bins, without #mu#mu CR"
