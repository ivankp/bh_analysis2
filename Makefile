STD := -std=c++14
DF := $(STD) -Isrc
CF := $(STD) -Wall -O3 -fmax-errors=3 -flto -Isrc
LF := $(STD) -flto

NPROC := $(shell nproc --all)

# ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_CFLAGS := -Wno-deprecated-declarations -pthread -m64 -I$(shell root-config --incdir)
ROOT_LIBS   := $(shell root-config --libs)

FJ_DIR    := $(shell fastjet-config --prefix)
FJ_CFLAGS := -I$(FJ_DIR)/include
FJ_LIBS   := -L$(FJ_DIR)/lib -lfastjet

C_hist_Hjets_mtop := $(ROOT_CFLAGS) $(FJ_CFLAGS) -DNPROC=$(NPROC)
L_hist_Hjets_mtop := $(ROOT_LIBS) -lTreePlayer $(FJ_LIBS)

SRC := src
BIN := bin
BLD := .build

SRCS := $(shell find $(SRC) -type f -name '*.cc')
DEPS := $(patsubst $(SRC)%.cc,$(BLD)%.d,$(SRCS))

GREP_EXES := grep -rl '^ *int \+main *(' $(SRC)
EXES := $(patsubst $(SRC)%.cc,$(BIN)%,$(shell $(GREP_EXES)))

NODEPS := clean
.PHONY: all clean

all: $(EXES)

#Don't create dependencies when we're cleaning, for instance
ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
-include $(DEPS)
endif

$(BIN)/test_re_axes: $(BLD)/re_axes.o
$(BIN)/hist_Hjets_mtop: $(BLD)/re_axes.o
$(BIN)/test_binning: $(BLD)/re_axes.o

$(DEPS): $(BLD)/%.d: $(SRC)/%.cc | $(BLD)
	$(CXX) $(DF) -MM -MT '$(@:.d=.o)' $< -MF $@

$(BLD)/%.o: | $(BLD)
	$(CXX) $(CF) $(C_$*) -c $(filter %.cc,$^) -o $@

$(BIN)/%: $(BLD)/%.o | $(BIN)
	$(CXX) $(LF) $(filter %.o,$^) -o $@ $(L_$*)

$(BLD) $(BIN):
	mkdir $@

clean:
	@rm -rfv $(BLD) $(BIN)
