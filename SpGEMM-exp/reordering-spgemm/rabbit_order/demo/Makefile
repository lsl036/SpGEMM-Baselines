GPERFTOOLS_DIR = /global/homes/r/raqib/rabbit_order/gperftools-2.9.1

CXXFLAGS += -Wall -Wextra -Wcast-align -Wcast-qual -Wconversion -Wfloat-equal \
	    -Wformat=2 -Winit-self -Wmissing-declarations \
	    -Wmissing-include-dirs -Wpointer-arith -Wredundant-decls \
	    -Wswitch-default -Wuninitialized -Wwrite-strings \
	    -Wno-sign-conversion -Wno-unused-function \
            -Wno-missing-declarations \
            -fopenmp -std=c++14 -mcx16 -O3 -DNDEBUG

CXXFLAGS += -I$(GPERFTOOLS_DIR)/include
LDFLAGS  += -L$(GPERFTOOLS_DIR)/lib
LDLIBS   += -ltcmalloc_minimal -lnuma
#LDLIBS   += -lnuma
TARGETS   = reorder

#SAVE_TO_FILE_FLAG = -DSAVE_TO_FILE
#
#ifeq ($(SAVE_TO_FILE), 1)
#	CXXFLAGS += $(SAVE_TO_FILE_FLAG)
#endif

$(TARGETS): %: %.cc
	$(LINK.cc) -MD -o $@ $< $(LDFLAGS) $(LDLIBS)
	@cp $*.d .$*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	      -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> .$*.P; \
	  rm -f $*.d

.PHONY: clean
clean:
	$(RM) $(TARGETS) $(TARGETS:%=.%.P)

-include .*.P

