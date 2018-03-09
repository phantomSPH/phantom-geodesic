#----------------------------------
#
#  Makefile for simple n-body code
#  Daniel Price, May 2016
#
#----------------------------------
FC=gfortran

### Old flags:
# FFLAGS=-O3 -Wall -Wextra -g -pedantic -std=f2008 -fcheck=all -Wno-compare-reals -fdefault-real-8 -fdefault-double-8#-freal-4-real-16#

### Phantom flags (+ a couple of extra):
FFLAGS= -O3 -Wall -Wextra -Wno-compare-reals -Wno-unused-dummy-argument -frecord-marker=4 -gdwarf-2 \
        -finline-functions-called-once -finline-limit=1500 -funroll-loops -ftree-vectorize \
        -std=f2008 -fall-intrinsics -fdefault-real-8 -fdefault-double-8
DEBUGFLAG= -g -fcheck=all -ffpe-trap=invalid,zero,overflow -finit-real=nan -finit-integer=nan -fbacktrace
ifeq ($(DEBUG), yes)
    FFLAGS += ${DEBUGFLAG}
    FFLAGS := $(FFLAGS:-O3=-O0)
endif

VPATH= ../src

# Default metric type
METRIC=schwarzschild
METRIC_FILE=metric_$(METRIC).f90

SRC=prompting.f90 inverse4x4.f90 $(METRIC_FILE) metric_tools.f90 io.f90 utils_testsuite.f90 utils_gr.f90 eos.f90 \
force_gr.f90 setup.f90 cons2prim.f90 step.f90 output.f90 test_metric.f90 checks.f90 test_cons2prim.f90\
 main.f90

ifeq ($(NEWTONIAN), yes)
	SRC=inverse4x4.f90 metric_minkowski.f90 metric_tools.f90 io.f90 utils_testsuite.f90 utils_gr.f90 eos.f90 \
	force_gr_newtonian.f90 setup.f90 cons2prim_newtonian.f90 step.f90 output.f90 test_metric.f90 checks.f90 test_cons2prim.f90 \
	main.f90
endif

OBJ=${SRC:.f90=.o}

SRC_TEST=${SRC:main.f90=test.f90}
OBJ_TEST=${SRC_TEST:.f90=.o}

%.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $@

grtest: ${OBJ}
	$(FC) $(FFLAGS) -o $@ $(OBJ)

test: $(OBJ_TEST)
	$(FC) $(FFLAGS) -o $@ $(OBJ_TEST)

clean:
	rm -f *.o *.mod
