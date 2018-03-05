include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG} ${FOMP}
   CFLAGS = -g ${FOMP}
else
   FFLAGS=${FOPT} ${FOMP}
   CFLAGS = -O3 ${FOMP}
endif

EXEC = neo

ifdef NEURAL_ROOT
   # neural net libraries
   NN=neo_neural.fann
else
   #stub file
   NN=neo_neural
endif

OBJECTS = neo_globals.o \
          neo_energy_grid.o \
          neo_interface.o \
          neo_allocate_profile.o \
          neo_umfpack.o \
          neo_sparse_solve.o \
          neo_equilibrium.o \
          neo_g_velocitygrids.o \
          neo_rotation.o \
          neo_nclass_dr.o \
          neo_theory.o \
          neo_transport.o \
          neo_3d_driver.o \
          neo_check.o \
          neo_compute_fcoll.o \
	  $(NN).o \
          neo_do.o \
          neo_error.o \
          neo_init.o \
          neo_init_serial.o \
          neo_make_profiles.o \
          neo_read_input.o \
          neo_run.o \
          neo_spitzer.o \
          st_to_cc.o

#-------------------------------------------------
# NEO_SOLVER can be set in platform Makefiles to:
#     slu   : SuperLU sparse solver library
#     petsc : PETSc library suite
#-------------------------------------------------
ifeq ($(NEO_SOLVER),slu)
   OBJECTS += neo_superlu.o
else
   OBJECTS += neo_superlu_stub.o
endif

.SUFFIXES : .o .f90 .F90

all: neo_lib.a $(EXEC).o $(EXTRA_LIBS)
	$(FC) $(FFLAGS) -o $(EXEC) $(EXEC).o neo_lib.a $(EXTRA_LIBS) $(LMATH) $(NN_LIB)

neo_lib.a: $(OBJECTS)
	$(ARCH) neo_lib.a $(OBJECTS)

.f90.o .F90.o :
	$(FC) $(FMATH) $(FFLAGS) $(NN_LIB) -c $< 

.c.o:
	$(CC) $(CFLAGS) $(CDEFS) $(EXTRA_HEADERS) -c $< $(VERBOSE)

clean:  
	rm -rf *.o neo_lib.a $(EXEC) *~
	cd ${GACODE_ROOT}/modules ; rm -f neo*.mod

-include ../install/make.ext.${GACODE_PLATFORM}
