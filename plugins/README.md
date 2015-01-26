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
	
	gen.sh                 for Gen level fitter.   
	
	accXreco.sh            for Efficiency calculation and fit.

	reco.sh                for Reco level fitter.
	


