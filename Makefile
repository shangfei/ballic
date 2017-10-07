#
# Files needed to compile the single component version of ballic
#
single_objects = ballic.single.o model.single.o grid.o

#
# Tipsy API
#
tipsy_objects = tipsy.o 

#
# Tillotson library
#
till_objects = tillotson/tillotson.o tillotson/tillinitlookup.o tillotson/tillsplint.o tillotson/interpol/coeff.o tillotson/interpol/interpol.o tillotson/interpol/brent.o tillotson/nr/nrcubicspline.o tillotson/nr/nrutil.o

#
# Icosahedron
#
fortran_objects = icosahedron.o
FC := gfortran

exe = ballic ballic.single ballic.multi modelsolve

#
# Code defines
#
#defs = -DTILL_PRESS_NP -DTILL_OUTPUT_ALL_WARNINGS -DTILL_PRESS_MELOSH

CFLAGS ?= -O3 -march=native 
FFLAGS ?= $(CFLAGS)

objects = $(single_objects)

default:
	@echo "Please specify, which tool to make (e.g., ballic, single, multi)."

all:
	default

# Standard version of ballic that reads an equilibrium model from a file and generates an particle distribution,
# that matches the desired density profile. Currently not reliably working for multi component models.
ballic: ballic.o $(tipsy_objects) $(fortran_objects)
	cc -o ballic ballic.o $(tipsy_objects) $(fortran_objects) -lm

# Single component version of ballic. It first calculates an equilibrium model for a given material and then
# generates a particle representation.
single: $(single_objects) $(tipsy_objects) $(till_objects) $(fortran_objects)
	cc -o ballic.single $(single_objects) $(tipsy_objects) $(till_objects) $(fortran_objects) -lm

# Two component version of ballic. It first calculates an equilibrium model for a given material and then
# generates a particle representation.
multi: ballic.multi.o $(tipsy_objects) $(till_objects) $(fortran_objects)
	cc -o ballic.multi ballic.multi.o $(tipsy_objects) $(till_objects) $(fortran_objects) -lm

# A two component model, that has a low density atmosphere of the mantle material
multi.atm: ballic.multi.atm.o $(tipsy_objects) $(till_objects) $(fortran_objects)
	cc -o ballic.multi.atm ballic.multi.atm.o $(tipsy_objects) $(till_objects) $(fortran_objects) -lm

#
# Calculates an equilibrium models for different rhoc and uc.
#
solve_model_array: solve_model_array.o $(till_objects)
	cc -o solve_model_array solve_model_array.o $(till_objects) -lm

#
# Calculates an equilibrium models for different rhoc and uc.
#
solve_model_array_multi: solve_model_array_multi.o $(till_objects)
	cc -o solve_model_array_multi solve_model_array_multi.o $(till_objects) -lm

#
# Calculates an equilibrium models for different rhoc and uc (parallel version).
#
solve_model_array_multi_omp: solve_model_array_multi.o $(till_objects)
	cc -o solve_model_array_multi_omp solve_model_array_multi.o $(till_objects) -lm -fopenmp

# Calculates equilibrium models for a given material and different initial densities and internal energies.
modelsolve: modelsolve.o $(till_objects) nr/gaussj.o
	cc -o modelsolve modelsolve.o $(till_objects) nr/gaussj.o -lm

#
# Calculates equilibrium models for a given material and different initial densities and internal energies
# using the structure equations as proposed in Alibert (2014).
#
modelsolve.alibert: modelsolve.alibert.o $(till_objects)
	cc -o modelsolve.alibert modelsolve.alibert.o $(till_objects) -lm

#
# Calculates equilibrium models for an ideal gas.
#
modelsolve.polytrope: modelsolve.polytrope.o $(till_objects)
	cc -o modelsolve.polytrope modelsolve.polytrope.o $(till_objects) -lm

#
# Solves the Lane-Emden equation.
#
lane-emden: lane-emden.o 
	cc -o lane-emden lane-emden.o -lm

clean:
	rm -f $(exe) $(objects)

cleanall:
	rm -f $(exe) $(objects) $(tipsy_objects) $(till_objects) $(fortran_objects)
