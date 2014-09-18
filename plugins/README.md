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



