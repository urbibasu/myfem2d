# -*- Makefile -*-
#
# Makefile for DynEarthSol3D
#
# Author: Eh Tan <tan2@earth.sinica.edu.tw>
#

## Execute "make" if making production run. Or "make opt=0 openmp=0" for debugging run.
##
## ndims = 3: 3D code; 2: 2D code
## opt = 1 ~ 3: optimized build; others: debugging build
## openmp = 1: enable OpenMP

ndims = 2
opt = 2

## Select C++ compiler
CXX = g++-mp-4.7

## Boost location and library name
BOOST_ROOT_DIR = /Users/eunseo/opt/boost_1_51_0

########################################################################
## Select compiler and linker flags
## (Usually you won't need to modify anything below)
########################################################################

BOOST_LDFLAGS = 
ifdef BOOST_ROOT_DIR
	BOOST_CXXFLAGS = -I$(BOOST_ROOT_DIR)
	BOOST_LDFLAGS += $(BOOST_ROOT_DIR)/stage/lib/libboost_program_options.a
endif

ifneq (, $(findstring g++, $(CXX))) # if using any version of g++
	CXXFLAGS = -g -std=c++0x
	LDFLAGS = -lm

	ifeq ($(opt), 1)
		CXXFLAGS += -O1
	else ifeq ($(opt), 2)
		CXXFLAGS += -O2
	else ifeq ($(opt), 3) # experimental, use at your own risk :)
		CXXFLAGS += -march=native -O3 -ffast-math -funroll-loops
	else # debugging flags
		CXXFLAGS += -O0 -Wall -Wno-unused-variable -Wno-unused-function -Wno-unknown-pragmas -fbounds-check -ftrapv
	endif

else ifneq (, $(findstring icpc, $(CXX))) # if using intel compiler, tested with v14
        CXXFLAGS = -g -std=c++0x
        LDFLAGS = -lm

        ifeq ($(opt), 1)
                CXXFLAGS += -O1
        else ifeq ($(opt), 2)
                CXXFLAGS += -O2
        else ifeq ($(opt), 3) # experimental, use at your own risk :)
                CXXFLAGS += -fast -fast-transcendentals -fp-model fast=2
        else # debugging flags
                CXXFLAGS += -O0 -check=uninit -check-pointers=rw -check-pointers-dangling=all -fp-trap-all=all
        endif

        ifeq ($(openmp), 1)
                CXXFLAGS += -fopenmp -DUSE_OMP
                LDFLAGS += -fopenmp
        endif

else
# the only way to display the error message in Makefile ...
all:
	@echo "Unknown compiler, check the definition of 'CXX' in the Makefile."
	@false
endif

##

SRCS =	\
	binaryio.cxx \
	main.cxx \
	fields.cxx \
	geometry.cxx \
	ic.cxx \
	input.cxx \
	matprops.cxx \
	mesh.cxx \
	output.cxx \
	rheology.cxx

INCS =	\
	array2d.hpp \
	binaryio.hpp \
	constants.hpp \
	parameters.hpp \
	matprops.hpp \
	utils.hpp \
	mesh.hpp \
	output.hpp

OBJS = $(SRCS:.cxx=.$(ndims)d.o)

EXE = myfem$(ndims)d


## Libraries

TRI_SRCS = triangle/triangle.c
TRI_INCS = triangle/triangle.h
TRI_OBJS = $(TRI_SRCS:.c=.o)

M_SRCS = $(TRI_SRCS)
M_INCS = $(TRI_INCS)
M_OBJS = $(TRI_OBJS)

## Action

.PHONY: all clean

all: $(EXE)

$(EXE): $(M_OBJS) $(OBJS)
	$(CXX) $(M_OBJS) $(OBJS) $(BOOST_LDFLAGS) $(LDFLAGS) -o $@

$(OBJS): %.$(ndims)d.o : %.cxx $(INCS)
	$(CXX) $(CXXFLAGS) $(BOOST_CXXFLAGS) -c $< -o $@

$(TRI_OBJS): %.o : %.c $(TRI_INCS)
	@# Triangle cannot be compiled with -O2
	$(CXX) $(CXXFLAGS) -O1 -DTRILIBRARY -DREDUCED -DANSI_DECLARATORS -c $< -o $@

deepclean:
	@rm -f $(TRI_OBJS) $(OBJS) $(EXE)

clean:
	@rm -f $(OBJS) $(EXE)

