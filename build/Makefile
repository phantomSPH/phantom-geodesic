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
FFLAGS= -O3 -Wall -Wextra -Wno-compare-reals -Wno-unused-dummy-argument -frecord-marker=4 \
        -finline-functions-called-once -finline-limit=1500 -funroll-loops -ftree-vectorize \
        -std=f2008 -fall-intrinsics -fdefault-real-8 -fdefault-double-8
DEBUGFLAG= -g -gdwarf-2 -fcheck=all -ffpe-trap=invalid,zero,overflow -finit-real=nan -finit-integer=nan -fbacktrace
OMPFLAGS = -fopenmp
OPENMP=no
ifeq ($(OPENMP), yes)
    FFLAGS += ${OMPFLAGS}
endif
ifeq ($(DEBUG), yes)
    FFLAGS += ${DEBUGFLAG}
    FFLAGS := $(FFLAGS:-O3=-O0)
endif

OMP_SCHEDULE=dynamic
OMP_NUM_THREADS=8
OMP_STACKSIZE=1024m

VPATH= ../src

# Default metric type
METRIC=kerr

ifeq ($(METRIC),bh)
    METRIC_FILE=metric_schwarzschild.f90
else ifeq ($(METRIC),flat)
    METRIC_FILE=metric_minkowski.f90
else
    METRIC_FILE=metric_$(METRIC).f90
endif

SETUP=geodesic

ifeq ($(SETUP), geodesic)
    SETUP_FILES=set_geodesic.f90 setup_geodesic.f90
else ifeq ($(SETUP), multigeodesics)
    SETUP_FILES=set_geodesic.f90 setup_multigeodesics.f90
else ifeq ($(SETUP), tde)
    SETUP_FILES=set_sphere.f90 setup_tde.f90
else ifeq ($(SETUP), multispheres)
    SETUP_FILES=set_sphere.f90 setup_multispheres.f90
else
    SETUP_FILES=setup_$(SETUP).f90
endif

SRC=options.f90 prompting.f90 utils.f90 utils_infiles.f90 inverse4x4.f90 $(METRIC_FILE) metric_tools.f90 io.f90 \
    utils_testsuite.f90 utils_gr.f90 eos.f90 force_gr.f90 $(SETUP_FILES) cons2prim.f90 energies.f90 step.f90 output.f90 \
    test_metric.f90 checks.f90 test_cons2prim.f90 init.f90 infile.f90 main.f90

ifeq ($(NEWTONIAN), yes)
	SRC=options.f90 prompting.f90 utils.f90 utils_infiles.f90 inverse4x4.f90 metric_minkowski.f90 metric_tools.f90 io.f90 \
   utils_testsuite.f90 utils_gr.f90 eos.f90 force_gr_newtonian.f90 $(SETUP_FILES) cons2prim_newtonian.f90 energies.f90 step.f90 output.f90 \
   test_metric.f90 checks.f90 test_cons2prim.f90 init.f90 infile.f90 main.f90
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
