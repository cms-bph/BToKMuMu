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

Run Fitter: ( more details in README.txt )
	
	mkdir plots/
	mkdir RootFiles/
	mkdir fitParameters/
	
	./fit <function> <input.root> <iBin>

Fit with batch mode using python: ( First of all, make sure that "test = False" in fit.py !!! )   
	
	./fit.py <function> <datatype> <input.root>  <iBin>  -b
	
	./fit.py angular_gen mc.gen ../RootFiles/MC_GENOnly/BToKMuMu_GENOnly_8TeV_genonly_v3-1.root  1  -b

The commands are collected in .sh files:  
	
	gen.sh                   for GEN level fitter.   
	
	accXreco.sh              for Efficiency calculation and fit.

	gen_R.sh                 for (RECO / Acceptance) angular fitting.

	reco.sh                  for RECO level fitter.

	angular2D_1a_Sm_local.sh for Bmass shape from signal MC, localy run.

	angular2D_1a_Sm.sh       for Bmass shape from signal MC, run with jobs.

	1b2a_jpsi.sh             for JpsiK peaking bkg.

	1b2a_psi.sh              for PsiK peaking bkg.

	prior.sh                 for Comb. bkg. from data sideband.
	
	test_d.sh                for data fit test with initial values.

	run.sh                   for data fit scanning with random initial values.

	fcn.sh                   for FCN values ploting and get the smallest one.

	ReFit.sh                 for final fitting results with best initial values.
	


