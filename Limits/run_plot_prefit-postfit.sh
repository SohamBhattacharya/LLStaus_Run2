#!/bin/bash

./plot_prefit-postfit_DisTauSF.py \
--input ~/nfs_dust/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/tmp/test_DisTauSF_mass_nbins1-60-80/DisTauSF/channels_all/eras_all/ZMT_wp-p*/postfit_s.root \
--outdir results/test_DisTauSF_mass_nbins1-60-80/channels_all/yields \
--era 2018 \
--channels mumu mutau

./plot_prefit-postfit_DisTauSF.py \
--input ~/nfs_dust/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/tmp/test_DisTauSF_mass_nbins2-60-70-80/DisTauSF/channels_all/eras_all/ZMT_wp-p*/postfit_s.root \
--outdir results/test_DisTauSF_mass_nbins2-60-70-80/channels_all/yields \
--era 2018 \
--channels mumu mutau

./plot_prefit-postfit_DisTauSF.py \
--input ~/nfs_dust/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/tmp/test_DisTauSF_mass_nbins1-60-80_no-dy-rateparam/DisTauSF/channels_all/eras_all/ZMT_wp-p*/postfit_s.root \
--outdir results/test_DisTauSF_mass_nbins1-60-80_no-dy-rateparam/channels_all/yields \
--era 2018 \
--channels mutau

./plot_prefit-postfit_DisTauSF.py \
--input ~/nfs_dust/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/tmp/test_DisTauSF_mass_nbins2-60-70-80_no-dy-rateparam/DisTauSF/channels_all/eras_all/ZMT_wp-p*/postfit_s.root \
--outdir results/test_DisTauSF_mass_nbins2-60-70-80_no-dy-rateparam/channels_all/yields \
--era 2018 \
--channels mutau
