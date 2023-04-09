
#inlcude make.inc of QE
include ../../make.inc

# for intel compiler
FFLAGS += -qopenmp -cpp
LDFLAGS += -qopenmp

# for gfortran compiler
#FFLAGS += -fopenmp -x f95-cpp-input 
#LDFLAGS += -fopenmp

#important: by default, Perturbo is compiled on top of QE-6.5
# uncomment the following line if you wanna compile it on top of QE-6.4.x
#FFLAGS += -D__QE64

#path to HDF5 library
IFLAGS += -I/projects/sg/brianr5/libraries/hdf5/include
HDF5_LIB = -L/projects/sg/brianr5/libraries/hdf5/lib -lhdf5 -lhdf5_fortran


# for debug mode
#FFLAGS += -fbacktrace -fbounds-check -fcheck=all 
#LDFLAGS += -fbacktrace -fbounds-check -fcheck=all 


MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../dft-d3 \
          $(MOD_FLAG)../../LR_Modules\
			 $(MOD_FLAG)../../PHonon/PH

          
PHMODS = ../../PHonon/PH/libph.a
LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../KS_Solvers/libks_solvers.a \
         ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a \
			../../UtilXlib/libutil.a ../../dft-d3/libdftd3qe.a

F90FLAGS = $(FFLAGS) $(FDFLAGS) $(MODFLAGS) $(IFLAGS)
LDFLAGS += $(HDF5_LIB)
