STD := -std=c++14
DF := $(STD) -Isrc
CF := $(STD) -Wall -Isrc
LF := $(STD)

ifneq (,${PREFIX})
PREFIX := ${PREFIX}
else
PREFIX := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
endif

ifeq (0, $(words $(findstring $(MAKECMDGOALS), rel)))
# development mode
CF += -O2 -g -fmax-errors=3
else
# release mode
CF += -O3 -flto -funroll-loops -march=native -mfpmath=sse
LF += -flto
endif

NPROC := $(shell nproc --all)

# ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_CFLAGS := -Wno-deprecated-declarations -pthread -m64
ROOT_CFLAGS += -I$(shell root-config --incdir)
ROOT_LIBS   := $(shell root-config --libs)

FJ_DIR    := $(shell fastjet-config --prefix)
FJ_CFLAGS := -I$(FJ_DIR)/include
FJ_LIBS   := -L$(FJ_DIR)/lib -lfastjet

LHAPDF_CFLAGS := $(shell lhapdf-config --cppflags)
LHAPDF_LIBS   := $(shell lhapdf-config --ldflags)

C_hist_Hjets := $(ROOT_CFLAGS) $(FJ_CFLAGS)
L_hist_Hjets := $(ROOT_LIBS) -lTreePlayer $(FJ_LIBS)

C_reweigh := $(ROOT_CFLAGS) $(LHAPDF_CFLAGS)
L_reweigh := $(ROOT_LIBS) $(LHAPDF_LIBS)

C_reweigh1 := $(ROOT_CFLAGS) $(LHAPDF_CFLAGS)
L_reweigh1 := $(ROOT_LIBS) $(LHAPDF_LIBS)

C_hist_Hjets_mtop := $(C_hist_Hjets)
L_hist_Hjets_mtop := $(L_hist_Hjets)

C_hist_example := $(C_hist_Hjets)
L_hist_example := $(L_hist_Hjets)

C_njets_test := $(C_hist_Hjets)
L_njets_test := $(L_hist_Hjets)

C_xsec_Hjets := $(ROOT_CFLAGS) $(FJ_CFLAGS)
L_xsec_Hjets := $(ROOT_LIBS) -lTreePlayer $(FJ_LIBS)

SRC := src
BIN := bin
BLD := .build

SRCS := $(shell find $(SRC) -type f -name '*.cc')
DEPS := $(patsubst $(SRC)%.cc,$(BLD)%.d,$(SRCS))

GREP_EXES := grep -rl '^ *int \+main *(' $(SRC)
EXES := $(patsubst $(SRC)%.cc,$(BIN)%,$(filter %.cc,$(shell $(GREP_EXES))))

HISTS := $(filter $(BIN)/hist_%,$(EXES))

NODEPS := clean
.PHONY: all clean

all: $(EXES)
rel: $(EXES)

#Don't create dependencies when we're cleaning, for instance
ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
-include $(DEPS)
endif

$(HISTS): $(BLD)/re_axes.o

$(DEPS): $(BLD)/%.d: $(SRC)/%.cc | $(BLD)
	$(CXX) $(DF) -MM -MT '$(@:.d=.o)' $< -MF $@

$(BLD)/%.o: | $(BLD)
	$(CXX) $(CF) $(C_$*) -DPREFIX="$(PREFIX)" -c $(filter %.cc,$^) -o $@

$(EXES): $(BIN)/%: $(BLD)/%.o | $(BIN)
	$(CXX) $(LF) $(filter %.o,$^) -o $@ $(L_$*)

$(BLD) $(BIN):
	mkdir $@

clean:
	@rm -rfv $(BLD) $(BIN)
