BAMTOOLS_ROOT=/srv/gsfs0/projects/curtis/ruping/tools/bamtools/
ZLIB_ROOT=/srv/gsfs0/projects/curtis/ruping/tools/zlib/current/
BOOST_ROOT=/srv/gsfs0/projects/curtis/ruping/tools/boost/current/
CXX=g++
BAMFLAGS=-lbamtools
CXXFLAGS=-lz
LBFLAGS=-Wl,-rpath,$(BAMTOOLS_ROOT)/lib/lib/:$(BOOST_ROOT)/lib
BOOSTFLAGS=-lboost_regex
PREFIX=$(CURDIR)
SRC=$(CURDIR)/src
TOOLSB=$(CURDIR)/utils/
BIN=/bin/
SOURCE_STA=Rseq_bam_stats.cpp
SOURCE_MFV=mappingFlankingVariants.cpp
SOURCE_REC=novelSnvFilter_ACGT.cpp
SOURCE_GS=grep_starts.cpp
STA=Rseq_bam_stats
MFV=mappingFlankingVariants
REC=novelSnvFilter_ACGT
GS=grep_starts

all: Rseq_bam_stats mappingFlankingVariants novelSnvFilter_ACGT grep_starts perl_scripts R_scripts lutils

.PHONY: all

Rseq_bam_stats:
	@mkdir -p $(PREFIX)/$(BIN)
	@echo "* compiling" $(SOURCE_STA)
	@$(CXX) $(SRC)/$(SOURCE_STA) -o $(PREFIX)/$(BIN)/$(STA) $(BAMFLAGS) $(CXXFLAGS) $(LBFLAGS) $(BOOSTFLAGS) -I $(BAMTOOLS_ROOT)/include/ -I $(ZLIB_ROOT)/include/ -I $(BOOST_ROOT)/include/ -L $(BAMTOOLS_ROOT)/lib/ -L $(ZLIB_ROOT)/lib/ -L $(BOOST_ROOT)/lib/

mappingFlankingVariants:
	@echo "* compiling" $(SOURCE_MFV)
	@$(CXX) $(SRC)/$(SOURCE_MFV) -o $(PREFIX)/$(BIN)/$(MFV) $(BAMFLAGS) $(CXXFLAGS) $(LBFLAGS) $(BOOSTFLAGS) -I $(BAMTOOLS_ROOT)/include/ -I $(ZLIB_ROOT)/include/ -I $(BOOST_ROOT)/include/ -L $(BAMTOOLS_ROOT)/lib/ -L $(ZLIB_ROOT)/lib/ -L $(BOOST_ROOT)/lib/

novelSnvFilter_ACGT:
	@echo "* compiling" $(SOURCE_REC)
	@$(CXX) $(SRC)/$(SOURCE_REC) -o $(PREFIX)/$(BIN)/$(REC) $(BAMFLAGS) $(CXXFLAGS) $(LBFLAGS) $(BOOSTFLAGS) -I $(BAMTOOLS_ROOT)/include/ -I $(ZLIB_ROOT)/include/ -I $(BOOST_ROOT)/include/ -L $(BAMTOOLS_ROOT)/lib/ -L $(ZLIB_ROOT)/lib/ -L $(BOOST_ROOT)/lib/

grep_starts:
	@echo "* compiling" $(SOURCE_GS)
	@$(CXX) $(SRC)/$(SOURCE_GS) -o $(PREFIX)/$(BIN)/$(GS) $(BAMFLAGS) $(CXXFLAGS) $(LBFLAGS) $(BOOSTFLAGS) -I $(BAMTOOLS_ROOT)/include/ -I $(ZLIB_ROOT)/include/ -I $(BOOST_ROOT)/include/ -L $(BAMTOOLS_ROOT)/lib/ -L $(ZLIB_ROOT)/lib/ -L $(BOOST_ROOT)/lib/

perl_scripts:
	@echo "* copying perl scripts"
	@cp $(SRC)/*.pl $(PREFIX)/$(BIN)/
	@echo "* done."

R_scripts:
	@echo "* copying R scripts"
	@cp $(SRC)/*.R $(PREFIX)/$(BIN)/
	@echo "* done."

lutils:
	@echo "* link utils"
	@ln -s $(TOOLSB)/* $(PREFIX)/$(BIN)/
	@echo "* done."

DTrace:
	@echo "* copying DTrace.pl"
	@cp $(SRC)/DTrace.pl $(PREFIX)/$(BIN)/
	@echo "* done."

clean:
	@echo "Cleaning up everthing."
	@rm -rf $(PREFIX)/$(BIN)/


.PHONY: clean
