#
# This is an example GNUmakefile for my packages
#
PACKAGE_NAME = Supera

# specific names for this package
SOURCES = $(wildcard *.cxx)
FMWK_HEADERS = LinkDef.h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))
HEADERS_DEST = larcv/app/$(PACKAGE_NAME)

# include options for this package
INCFLAGS  = -I.                       #Include itself
INCFLAGS += $(shell larcv-config --includes)

ifdef LARCV_EDEPSIM_INCLUDES
INCFLAGS += $(LARCV_EDEPSIM_INCLUDES)
LDFLAGS  += $(LARCV_EDEPSIM_LIBS)
endif

#MINE ??????????
#INCFLAGS += -I/app/edep/include
#INCFLAGS += -I/app/geant4/include/Geant4
#

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

include $(LARCV_BASEDIR)/Makefile/Makefile.${OSNAME}

LDFLAGS += $(shell larcv-config --libs)

# call the common GNUmakefile
include $(LARCV_BASEDIR)/Makefile/GNUmakefile.CORE



pkg_build:
pkg_clean:
