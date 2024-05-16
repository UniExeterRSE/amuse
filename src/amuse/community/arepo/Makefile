# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include $(AMUSE_DIR)/config.mk

# File populating AREPOFLAGS for Arepo ifdefs
include Arepoflags.mk

HDF5_FLAGS += -DH5_USE_16_API

MPICXX   ?= mpicxx
CXX = $(MPICXX)
CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) $(GSL_FLAGS) $(AREPOFLAGS) $(HDF5_FLAGS)
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = src/libarepo.a

all: arepo_worker 

clean:
	$(RM) -rf __pycache__
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ arepo_worker worker_code.cc
	make -C src clean

distclean: clean
	make -C src distclean

$(CODELIB): .FORCE
	export AREPOFLAGS
	make -C src all

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py ArepoInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py ArepoInterface -o $@

arepo_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(GSL_FLAGS) $(HDF5_FLAGS) $< $(OBJS) $(CODELIB) -o $@ $(GMP_LIBS) $(GSL_LIBS) $(HDF5_LIBS)

.cc.o: $<
	$(MPICXX) $(CXXFLAGS) $(CODELIB) -c -o $@ $< 

.FORCE:
