#----------------------------------
#
#  Makefile for simple n-body code
#  Daniel Price, May 2016
#
#----------------------------------
FC=gfortran
FFLAGS=-O3 -Wall -Wextra -pedantic -std=f2008 -fdefault-real-8 -fdefault-double-8

# Default metric type
METRIC=schwarzschild
METRIC_FILE=metric_$(METRIC).f90

SRC=$(METRIC_FILE) io.f90 utils_testsuite.f90 utils_gr.f90 eos.f90 \
init.f90 cons2prim.f90 force_gr.f90 step.f90 output.f90 test_metric.f90 checks.f90 test_cons2prim.f90\
 main.f90
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
