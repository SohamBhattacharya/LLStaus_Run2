#ERA=2016_preVFP
#python -u -m pepper.runproc \
#processor_MET-trigger-efficiency.py configs/config_trigger_eff_MET_$ERA.json \
#-o ./output/pepper_condor/MET-trigger-efficiency/$ERA \
#--statedata ./output/pepper_condor/MET-trigger-efficiency/$ERA/pepper_condor_fake_rate.coffea \
#-i /nfs/dust/cms/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Analysis/configs/setup_env_sobhatta.sh \
#--condor 2000 \
#--retries 10


ERA=2016_postVFP
python -u -m pepper.runproc \
processor_MET-trigger-efficiency.py configs/config_trigger_eff_MET_$ERA.json \
-o ./output/pepper_condor/MET-trigger-efficiency/$ERA \
--statedata ./output/pepper_condor/MET-trigger-efficiency/$ERA/pepper_condor_fake_rate.coffea \
-i /nfs/dust/cms/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Analysis/configs/setup_env_sobhatta.sh \
--condor 2000 \
--retries 10


ERA=2018
#python -u -m pepper.runproc \
#processor_MET-trigger-efficiency.py configs/config_trigger_eff_MET_$ERA.json \
#-o ./output/pepper_condor/MET-trigger-efficiency/$ERA \
#--statedata ./output/pepper_condor/MET-trigger-efficiency/$ERA/pepper_condor_fake_rate.coffea \
#-i /nfs/dust/cms/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Analysis/configs/setup_env_sobhatta.sh \
#--condor 2000 \
#--retries 10