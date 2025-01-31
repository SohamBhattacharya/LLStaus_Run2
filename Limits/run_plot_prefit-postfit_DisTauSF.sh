#!/bin/bash -x

ERA="$1"

#./plot_prefit-postfit_DisTauSF.py \
#--input ~/nfs_dust/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/tmp/test_DisTauSF_mass_nbins1-60-80/DisTauSF/channels_all/eras_all/ZMT_wp-p*/postfit_s.root \
#--outdir results/test_DisTauSF_mass_nbins1-60-80/channels_all/yields \
#--era ${ERA} \
#--channels mumu mutau
#
#./plot_prefit-postfit_DisTauSF.py \
#--input ~/nfs_dust/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/tmp/test_DisTauSF_mass_nbins2-60-70-80/DisTauSF/channels_all/eras_all/ZMT_wp-p*/postfit_s.root \
#--outdir results/test_DisTauSF_mass_nbins2-60-70-80/channels_all/yields \
#--era ${ERA} \
#--channels mumu mutau

./plot_prefit-postfit_DisTauSF.py \
--input /home/soham/mnt/desy_dust/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/results/DisTauSF/DisTauSF_mass_nbins1-60-80/DisTauSF/channels_all/eras_${ERA}/ZMT_wp-p*/postfit_s.root \
--outdir results/DisTauSF/DisTauSF_mass_nbins1-60-80/DisTauSF/channels_all/eras_${ERA}/prefit-postfit \
--era ${ERA} \
--channels mutau

#./plot_prefit-postfit_DisTauSF.py \
#--input ~/nfs_dust/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Limits/tmp/test_DisTauSF_mass_nbins2-60-70-80_no-dy-rateparam/DisTauSF/channels_all/eras_all/ZMT_wp-p*/postfit_s.root \
#--outdir results/test_DisTauSF_mass_nbins2-60-70-80_no-dy-rateparam/channels_all/yields \
#--era ${ERA} \
#--channels mutau
