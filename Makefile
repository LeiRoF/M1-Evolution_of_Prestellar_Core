

# Choose compiler
 FC = gfortran

# Options and Path
 FFLAGS = -g -fno-second-underscore -Wall -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow
# FFLAGS = -O3


##############################################################################
#----------------------------------------------------------------------------

.MAKEOPTS: -k -s

.SUFFIXES: .f90 .f

.f90.o:
	$(FC) $(FFLAGS) -c $*.f90



OBJ    = Constants.o \
         Prestel.o

# List all the "rules" that must be executed when "make" is typed in the command line
ALL: Test PartialClean
	echo "!!!        Compilation OK        !!!"

# Rule to compile the main and all the modules
Test : $(OBJ) main_test.o
	$(FC) $(FFLAGS) main_test.o $(OBJ) -o main_test

# Rule to remove all the intermediate files
PartialClean:
	rm *.o *mod

# Rule to remove all the intermediate files and the executable
clean:
	rm *.o BGM *mod

##############################################################################
#----------------------------------------------------------------------------

