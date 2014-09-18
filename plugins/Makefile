GCCFLAGS=-Wall -g
ROOTFLAGS=`root-config --cflags --libs`

<<<<<<< HEAD
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
=======
SELECTORFLAGS=$(GCCFLAGS) $(ROOTFLAGS) -lProof
FIGURESFLAGS=$(GCCFLAGS) $(ROOTFLAGS)
FITFLAGS=$(GCCFLAGS) $(ROOTFLAGS) -lRooFit -lRooFitCore

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
>>>>>>> 079f0fc10915ff5ec232d76f56b47756c1255c32
