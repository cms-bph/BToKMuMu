##########################
# B To K Mu Mu Analysis 
##########################

##########
## Build 
##########

>  cd /path/to/your/destinaiton

>  cmsrel CMSSW_5_3_20

>  cd CMSSW_5_3_20/src/

>  cmsenv

>  cp -r /afs/cern.ch/user/g/gechen/gechen/public/BToKMuMu/RecoVertex ./

>  scram b

>  mkdir ana

>  cd ana/

>  git clone https://github.com/cms-bph/BToKMuMu.git 

>  cd ../

>  scram b

########
## Run 
########

>  cd ana/BToKMuMu/python

>  cmsRun btokmumu_cfg.py 


