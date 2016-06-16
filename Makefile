#----------------------------------
#
#  Makefile for simple n-body code
#  Daniel Price, May 2016
#
#----------------------------------
FC=gfortran
FFLAGS=-O3 -Wall -Wextra -pedantic -std=f2008 -fdefault-real-8 -fdefault-double-8

SRC=utils_gr.f90 init.f90 metric_schwarzschild.f90 cons2prim.f90 step.f90 output.f90 checks.f90 test.f90
OBJ=${SRC:.f90=.o}

%.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $@

grtest: ${OBJ}
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean:
	rm -f *.o *.mod
