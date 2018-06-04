SHELL = bash
CXX := g++
STD := -std=c++14
CPPFLAGS := $(STD) -Iinclude
CXXFLAGS := $(STD) -Wall -O3 -flto -Iinclude -fmax-errors=3
# CXXFLAGS := $(STD) -Wall -g -Iinclude -fmax-errors=3
LDFLAGS := $(STD) -O3 -flto
# LDFLAGS :=
LDLIBS :=

ifeq (,${PREFIX})
PREFIX := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
endif
CXXFLAGS += -DPREFIX="$(PREFIX)"

SRC := src
BIN := bin
BLD := .build
EXT := .cc

.PHONY: all clean

ifeq (0, $(words $(findstring $(MAKECMDGOALS), clean)))

LDFLAGS += $(shell sed -r 's/([^:]+)(:|$$)/ -L\1/g' <<< "$$LIBRARY_PATH")

ROOT_CXXFLAGS := $(shell root-config --cflags | sed 's/ -std=c++[^ ]\+ / /')
ROOT_LDLIBS   := $(shell root-config --libs)

FJ_PREFIX   := $(shell fastjet-config --prefix)
FJ_CXXFLAGS := -I$(FJ_PREFIX)/include
FJ_LDLIBS   := -L$(FJ_PREFIX)/lib -Wl,-rpath=$(FJ_PREFIX)/lib -lfastjet

LHAPDF_PREFIX   := $(shell lhapdf-config --prefix)
LHAPDF_CXXFLAGS := $(shell lhapdf-config --cppflags)
LHAPDF_LDLIBS   := $(shell lhapdf-config --libs) -Wl,-rpath=$(LHAPDF_PREFIX)/lib

# RPATH
rpath_script := ldd $(shell root-config --libdir)/libTreePlayer.so \
  | sed -nr 's|.*=> (.+)/.+\.so[.0-9]* \(.*|\1|p' \
  | sort -u \
  | sed -nr '/^(\/usr)?\/lib/!s/^/-Wl,-rpath=/p'
ROOT_LDLIBS += $(shell $(rpath_script))

C_check_tree := $(ROOT_CXXFLAGS)
L_check_tree := $(ROOT_LDLIBS) -lTreePlayer

C_reweigh := $(ROOT_CXXFLAGS)
L_reweigh := $(ROOT_LDLIBS) $(LHAPDF_LDLIBS) -lboost_program_options

C_uncert := $(ROOT_CXXFLAGS)
L_uncert := $(ROOT_LDLIBS) $(LHAPDF_LDLIBS)

C_reweighter := $(ROOT_CXXFLAGS) $(LHAPDF_CXXFLAGS)

C_reweigh1 := $(ROOT_CXXFLAGS) $(LHAPDF_CXXFLAGS)
L_reweigh1 := $(ROOT_LDLIBS) $(LHAPDF_LDLIBS)

C_reweigh_threaded := $(ROOT_CXXFLAGS) $(LHAPDF_CXXFLAGS)
L_reweigh_threaded := $(ROOT_LDLIBS) $(LHAPDF_LDLIBS)

C_dep_scale := $(ROOT_CXXFLAGS) $(FJ_CXXFLAGS)
L_dep_scale := $(ROOT_LDLIBS) $(LHAPDF_LDLIBS) $(FJ_LDLIBS)

C_dep_R_scale := $(ROOT_CXXFLAGS) $(FJ_CXXFLAGS)
L_dep_R_scale := $(ROOT_LDLIBS) $(LHAPDF_LDLIBS) $(FJ_LDLIBS)

C_hist_Hjets := $(ROOT_CXXFLAGS) $(FJ_CXXFLAGS)
L_hist_Hjets := $(ROOT_LDLIBS) -lTreePlayer $(FJ_LDLIBS)

C_hist_Hjets_yy := $(C_hist_Hjets)
L_hist_Hjets_yy := $(L_hist_Hjets)

C_hist_Hjets_mtop := $(C_hist_Hjets)
L_hist_Hjets_mtop := $(L_hist_Hjets)

C_hist_Hjets_isolation := $(C_hist_Hjets)
L_hist_Hjets_isolation := $(L_hist_Hjets)

C_hist_hgam := $(C_hist_Hjets)
L_hist_hgam := $(L_hist_Hjets)

C_hist_hgam_sb := $(C_hist_Hjets)
L_hist_hgam_sb := $(L_hist_Hjets)

C_hist_example := $(C_hist_Hjets)
L_hist_example := $(L_hist_Hjets)

C_njets_test := $(C_hist_Hjets)
L_njets_test := $(L_hist_Hjets)

C_xsec_Hjets := $(ROOT_CXXFLAGS) $(FJ_CXXFLAGS)
L_xsec_Hjets := $(ROOT_LDLIBS) -lTreePlayer $(FJ_LDLIBS)

C_Higgs2diphoton := $(ROOT_CXXFLAGS)
C_test_H2AA := $(ROOT_CXXFLAGS)
L_test_H2AA := $(ROOT_LDLIBS)

C_test_binner := $(ROOT_CXXFLAGS)
L_test_binner := $(ROOT_LDLIBS) -lTreePlayer

C_unweighted := $(ROOT_CXXFLAGS) $(FJ_CXXFLAGS)
L_unweighted := $(ROOT_LDLIBS) -lTreePlayer $(FJ_LDLIBS)

C_unweighted2text := $(ROOT_CXXFLAGS)
L_unweighted2text := $(ROOT_LDLIBS)

C_hist_Hjets_ang := $(C_hist_Hjets)
L_hist_Hjets_ang := $(L_hist_Hjets)

C_fit_Hjets_ang := $(ROOT_CXXFLAGS)
L_fit_Hjets_ang := $(ROOT_LDLIBS)

C_var_Hjets_angles := $(ROOT_CXXFLAGS) $(FJ_CXXFLAGS)
L_var_Hjets_angles := $(ROOT_LDLIBS) -lTreePlayer $(FJ_LDLIBS)

C_var_H1j_4mom := $(ROOT_CXXFLAGS)
L_var_H1j_4mom := $(ROOT_LDLIBS) -lTreePlayer

C_csv_4mom := $(ROOT_CXXFLAGS)
L_csv_4mom := $(ROOT_LDLIBS) -lTreePlayer -lboost_iostreams

C_merge_nlo := $(ROOT_CXXFLAGS)
L_merge_nlo := $(ROOT_LDLIBS)

C_hist_Hj_angular := $(ROOT_CXXFLAGS) $(FJ_CXXFLAGS) -DIMPL=hist/angular.hh
L_hist_Hj_angular := $(ROOT_LDLIBS) -lTreePlayer $(FJ_LDLIBS)

C_hist_Hj_hgam_sb := $(ROOT_CXXFLAGS) $(FJ_CXXFLAGS) -DIMPL=hist/hgam_sb.hh
L_hist_Hj_hgam_sb := $(ROOT_LDLIBS) -lTreePlayer $(FJ_LDLIBS)

SRCS := $(shell find $(SRC) -type f -name '*$(EXT)' \
  -not -name 'hist_Hj$(EXT)')
DEPS := $(SRCS:$(SRC)/%$(EXT)=$(BLD)/%.d)

GREP_EXES := grep -rl '^[[:blank:]]*int \+main *(' $(SRC) \
  --include='*$(EXT)' --exclude='hist_Hj$(EXT)'
EXES := $(patsubst $(SRC)%$(EXT),$(BIN)%,$(shell $(GREP_EXES)))

HISTS := $(filter $(BIN)/hist_%,$(EXES))

HISTS2_SUF := $(patsubst include/hist/%.hh,%,$(wildcard include/hist/*.hh))
DEPS2 := $(HISTS2_SUF:%=$(BLD)/hist_Hj_%.d)
HISTS2 := $(HISTS2_SUF:%=$(BIN)/hist_Hj_%)
EXES += $(HISTS2)

all: $(EXES)

$(HISTS): $(BLD)/re_axes.o

$(BIN)/reweigh $(BIN)/dep_scale $(BIN)/dep_R_scale \
: $(BLD)/reweighter.o

$(BIN)/test_H2AA $(BIN)/hist_Hjets_isolation \
$(BIN)/hist_hgam $(BIN)/hist_hgam_sb \
$(BIN)/hist_Hjets_yy $(BIN)/unweighted $(BIN)/hist_Hjets_ang \
$(BIN)/var_Hjets_angles \
$(BIN)/var_H1j_4mom $(BIN)/csv_4mom \
: $(BLD)/Higgs2diphoton.o

$(BIN)/unweighted $(BIN)/unweighted2text $(BIN)/hist_Hjets_ang \
$(BIN)/fit_Hjets_ang $(BIN)/var_Hjets_angles \
$(BIN)/var_H1j_4mom $(BIN)/csv_4mom \
: $(BLD)/program_options.o

$(HISTS2): $(BLD)/re_axes.o $(BLD)/Higgs2diphoton.o $(BLD)/program_options.o

-include $(DEPS) $(DEPS2)

.SECONDEXPANSION:

$(DEPS): $(BLD)/%.d: $(SRC)/%$(EXT) | $(BLD)/$$(dir %)
	$(CXX) $(CPPFLAGS) $(C_$*) -MM -MT '$(@:.d=.o)' $< -MF $@

$(BLD)/%.o: | $(BLD)
	$(CXX) $(CXXFLAGS) $(C_$*) -c $(filter %$(EXT),$^) -o $@

$(BIN)/%: $(BLD)/%.o | $(BIN)
	$(CXX) $(LDFLAGS) $(filter %.o,$^) -o $@ $(LDLIBS) $(L_$*)

$(DEPS2): $(BLD)/hist_Hj_%.d: $(SRC)/hist_Hj$(EXT) include/hist/%.hh | $(BLD)/$$(dir %)
	$(CXX) $(CPPFLAGS) $(C_$*) -DIMPL=hist/$*.hh -MM -MT '$(@:.d=.o)' $< -MF $@

$(BIN):
	mkdir -p $@

$(BLD)/%/:
	mkdir -p $@

endif

clean:
	@rm -rfv $(BLD) $(BIN)
