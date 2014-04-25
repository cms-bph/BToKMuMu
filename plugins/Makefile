GCCFLAGS=-Wall -g
ROOTFLAGS=`root-config --cflags --libs`

TOOLSFLAGS=$(GCCFLAGS) -I$(ROOTSYS)/include -c
SELECTORFLAGS=$(GCCFLAGS) $(ROOTFLAGS) -lProof
#FIGURESFLAGS=$(GCCFLAGS) $(ROOTFLAGS)
#FITFLAGS=$(GCCFLAGS) $(ROOTFLAGS) -lRooFit -lRooFitCore -lMinuit

all: sel
	@echo "Full build successful."

sel: KMUMU.cc
	g++ $(SELECTORFLAGS) $< -o $@


clean:
	rm -rf  sel 
