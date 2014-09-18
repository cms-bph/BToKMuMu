GCCFLAGS=-Wall -g
ROOTFLAGS=`root-config --cflags --libs`

TOOLSFLAGS=$(GCCFLAGS) -I$(ROOTSYS)/include -c
SELECTORFLAGS=$(GCCFLAGS) $(ROOTFLAGS) -lProof
FIGURESFLAGS=$(GCCFLAGS) $(ROOTFLAGS)
FITFLAGS=$(GCCFLAGS) $(ROOTFLAGS) -lRooFit -lRooFitCore -lMinuit

all: sel fig fit
	@echo "Full build successful."


sel: SingleBuToKMuMuSelector.cc
	g++ $(SELECTORFLAGS) $< -o $@

fig: SingleBuToKMuMuFigures.cc
	g++ $(FIGURESFLAGS) $< -o $@

fit: SingleBuToKMuMuFitters.cc
	g++ $(FITFLAGS) $< -o $@	

clean:
	rm -rf sel fig fit 
