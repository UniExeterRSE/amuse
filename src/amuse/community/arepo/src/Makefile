MPICXX ?= mpicxx
MPICC ?= mpicc
# TODO: Determine whether this is needed as included in arepo/Makefile.
ifeq ($(origin AMUSE_DIR), undefined)
	AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include $(AMUSE_DIR)/config.mk

CC       =  $(MPICC)  # sets the C-compiler

# GSL_FLAGS  = -I/opt/Homebrew/include  # TODO: Need to make GSL_FLAGS generalisable.
CFLAGS   += -Wall -g $(GSL_FLAGS)
CXXFLAGS += $(CFLAGS) 

HDF5_FLAGS += -DH5_USE_16_API

LDFLAGS  += -lm $(MUSE_LD_FLAGS) $(HDF5_FLAGS)
LIBS += $(HDF5_LIBS)

CODELIB = libarepo.a

SUBDIR = add_backgroundgrid cooling debug_md5 domain fof gitversion gravity \
	hydro init io mesh mesh/voronoi mpi_utils ngbtree star_formation subfind \
	time_integration utils
SRCS = $(foreach fd, $(SUBDIR), $(wildcard $(fd)/*.c))
SRCS += main/allvars.c

CODEOBJS = $(SRCS:c=o)
CODEOBJS += main/allvars.o main/run.o main/main.o

# Arepoflags.mk populates AREPOFLAGS for Arepo ifdefs.
#
# These are defined in this file instead of being added to arepoflags.h by
# config-makefile and prepare-config.perl from a Config.sh (which would be
# typical of arepo).
include ../Arepoflags.mk

CXXFLAGS += $(AREPOFLAGS)

# Add includes
CXXFLAGS += $(HDF5_FLAGS)

AR = ar ruv
RANLIB = ranlib
RM = rm

all: $(CODELIB)

clean:
	$(RM) -f *.o *.a
	$(RM) $(SRCS:c=o)

distclean: clean

$(CODELIB): $(CODEOBJS)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS)
	$(RANLIB) $@

.cc.o: $<
	$(MPICXX) $(CXXFLAGS) -c -o $@ $< $(LIBS) 

.c.o: $<
	$(MPICC) $(CXXFLAGS) -c -o $@ $< $(LIBS) 
