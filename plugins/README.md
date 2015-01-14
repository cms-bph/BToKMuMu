Instructions for Plugins (run2012v0)

lxplus6 environment:   
	
	export VO_CMS_SW_DIR=/opt/exp_soft/cms
	source /afs/cern.ch/sw/lcg/external/gcc/4.7/x86_64-slc6-gcc47-opt/setup.sh
	version=5.34.14
	source /afs/cern.ch/sw/lcg/app/releases/ROOT/${version}/x86_64-slc6-gcc47-opt/root/bin/thisroot.sh

Build:
	
	make

Run Selection:
	
	./sel <datatype> <cut> <input.root> <output.root>

Run Fitter:
	
	mkdir plots/
	mkdir RootFiles/
	mkdir fitParameters/
	
	./fit <function> <input.root>

0. Test on Data Bmass fitter            >>>>>>>>>>>>>>>

   source Fitter_bmass_Data_AllCut.sh

1. Angular fitter in gen level          >>>>>>>>>>>>>>>

   ./fit angular_gen BToKMuMu_GENOnly_8TeV_genonly_v3-1.root

2. Calculate efficiency                 >>>>>>>>>>>>>>>

   ./fit accXrecoEff AllCut/MC_Signal_AllCuts_v3.root

3. Angular fitter in reco level         >>>>>>>>>>>>>>>

   ./fit angular_reco AllCut/MC_Signal_AllCuts_v3.root

4. Get bmass shape from signal mc       >>>>>>>>>>>>>>>

   ./fit angular2D_1a_Sm AllCut/MC_Signal_AllCuts_v3.root

5. Get bmass peaking bkg. from mc       >>>>>>>>>>>>>>>

   ./fit angular2D_1b_YpPm_Jpsi AllCut/MC_JPsiK_AllCuts_v3.root
   ./fit angular2D_1b_YpPm_Psi AllCut/MC_Psi2SK_AllCuts_v3.root

6. Get cosThetaL peaking bkg. from mc   >>>>>>>>>>>>>>>

   ./fit angular2D_2a_PkPl_Jpsi AllCut/MC_JPsiK_AllCuts_v3.root
   ./fit angular2D_2a_PkPl_Psi AllCut/MC_Psi2SK_AllCuts_v3.root

7. Get cosThetaL comb. bkg. from data   >>>>>>>>>>>>>>>

   ./fit angular2D_prior AllCut/Data_KMuMu_AllCuts_v3.root

8. Angular fitter in data               >>>>>>>>>>>>>>>

   ./fit angular2D AllCut/Data_KMuMu_AllCuts_v3.root 1




