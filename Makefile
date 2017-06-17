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

# RPATH
rpath_script := ldd `root-config --libdir`/libTreePlayer.so \
  | sed -n 's/.*=> \(.*\)\/.\+\.so[^ ]* (.*/\1/p' \
  | sort | uniq \
  | sed '/^\/lib/d;/^\/usr\/lib/d' \
  | sed 's/^/-Wl,-rpath=/'
ROOT_LIBS += $(shell $(rpath_script))

FJ_DIR    := $(shell fastjet-config --prefix)
FJ_CFLAGS := -I$(FJ_DIR)/include
FJ_LIBS   := -L$(FJ_DIR)/lib -lfastjet -Wl,-rpath=$(FJ_DIR)/lib

LHAPDF_DIR    := $(shell fastjet-config --prefix)
LHAPDF_CFLAGS := $(shell lhapdf-config --cppflags)
LHAPDF_LIBS   := $(shell lhapdf-config --ldflags) -Wl,-rpath=$(LHAPDF_DIR)/lib

C_check_tree := $(ROOT_CFLAGS)
L_check_tree := $(ROOT_LIBS) -lTreePlayer

C_reweigh := $(ROOT_CFLAGS)
L_reweigh := $(ROOT_LIBS) $(LHAPDF_LIBS) -lboost_program_options

C_uncert := $(ROOT_CFLAGS)
L_uncert := $(ROOT_LIBS) $(LHAPDF_LIBS)

C_reweighter := $(ROOT_CFLAGS) $(LHAPDF_CFLAGS)

C_reweigh1 := $(ROOT_CFLAGS) $(LHAPDF_CFLAGS)
L_reweigh1 := $(ROOT_LIBS) $(LHAPDF_LIBS)

C_reweigh_threaded := $(ROOT_CFLAGS) $(LHAPDF_CFLAGS)
L_reweigh_threaded := $(ROOT_LIBS) $(LHAPDF_LIBS)

C_dep_scale := $(ROOT_CFLAGS) $(FJ_CFLAGS)
L_dep_scale := $(ROOT_LIBS) $(LHAPDF_LIBS) $(FJ_LIBS)

C_dep_R_scale := $(ROOT_CFLAGS) $(FJ_CFLAGS)
L_dep_R_scale := $(ROOT_LIBS) $(LHAPDF_LIBS) $(FJ_LIBS)

C_hist_Hjets := $(ROOT_CFLAGS) $(FJ_CFLAGS)
L_hist_Hjets := $(ROOT_LIBS) -lTreePlayer $(FJ_LIBS)

C_hist_Hjets_mtop := $(C_hist_Hjets)
L_hist_Hjets_mtop := $(L_hist_Hjets)

C_hist_Hjets_isolation := $(C_hist_Hjets)
L_hist_Hjets_isolation := $(L_hist_Hjets)

C_hist_hgam := $(C_hist_Hjets)
L_hist_hgam := $(L_hist_Hjets)

C_hist_example := $(C_hist_Hjets)
L_hist_example := $(L_hist_Hjets)

C_njets_test := $(C_hist_Hjets)
L_njets_test := $(L_hist_Hjets)

C_xsec_Hjets := $(ROOT_CFLAGS) $(FJ_CFLAGS)
L_xsec_Hjets := $(ROOT_LIBS) -lTreePlayer $(FJ_LIBS)

C_Higgs2diphoton := $(ROOT_CFLAGS)
C_test_H2AA := $(ROOT_CFLAGS)
L_test_H2AA := $(ROOT_LIBS)

C_test_binner := $(ROOT_CFLAGS)
L_test_binner := $(ROOT_LIBS) -lTreePlayer

SRC := src
BIN := bin
BLD := .build

SRCS := $(shell find $(SRC) -type f -name '*.cc')
DEPS := $(patsubst $(SRC)%.cc,$(BLD)%.d,$(SRCS))

GREP_EXES := grep -rl '^[[:blank:]]*int \+main *(' $(SRC)
EXES := $(patsubst $(SRC)%.cc,$(BIN)%,$(shell $(GREP_EXES)))

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
$(BIN)/reweigh $(BIN)/dep_scale $(BIN)/dep_R_scale: $(BLD)/reweighter.o
$(BIN)/test_H2AA $(BIN)/hist_Hjets_isolation $(BIN)/hist_hgam: \
  $(BLD)/Higgs2diphoton.o

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
