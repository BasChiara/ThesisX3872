ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs) -lTMVA -lRooFit -lRooFitCore -lMinuit

CXX           = g++
CXXFLAGS      = -g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2
LD            = g++
LDFLAGS       = -g -lGenVector # con errori tipo 'undefined reference to `ROOT::Math::GenVector::Throw(char const*)'
SOFLAGS       = -shared


ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)


CXXFLAGS      += $(ROOTCFLAGS)
#CXX           += -I./
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) 
GLIBS          = $(filter-out -lNew , $(NGLIBS))

INCLUDEDIR       = ./
CXX	         += -I$(INCLUDEDIR) -I.
OUTLIB	         = $(INCLUDEDIR)/lib/

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ./lib/


$(OUTLIB)Prep.o: $(INCLUDEDIR)src/PrepAnalysis.C 
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)PrepAnalysis.o $<
$(OUTLIB)CheckGenLev.o: $(INCLUDEDIR)src/CheckGenLev.C 
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)CheckGenLev.o $<
$(OUTLIB)CheckGenLevPsi2S.o: $(INCLUDEDIR)src/CheckGenLevPsi2S.C 
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)CheckGenLevPsi2S.o $<
# ==================== GEN LEVEL CHECK =============================================
GenLev:$(INCLUDEDIR)macro.cc\
		$(OUTLIB)Prep.o\
		$(OUTLIB)CheckGenLev.o\
		$(OUTLIB)CheckGenLevPsi2S.o
			$(CXX) $(CXXFLAGS) -ldl -o analysisGenLevel $(OUTLIB)/*.o $(GLIBS) $(LDFLAGS) $ $<
GenLev.clean:
			 rm -f GenLev 

clean:
		rm -f $(OUTLIB)*.o
		rm -f GenLev

gen: GenLev