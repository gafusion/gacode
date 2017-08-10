include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

neo_libs= -L$(GACODE_ROOT)/shared/harvest_client/ \
              -I$(GACODE_ROOT)/shared/harvest_client/ \
              -lharvest

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC = neo

ifdef FANN_ROOT
    # neural-net file
    JBSNN=neo_jbs_nn.fann
    JBSNN_LIBS=-I$(GACODE_ROOT)/shared/neural/ -L$(GACODE_ROOT)/shared/neural/ -lbrainfuse -L${FANN_ROOT}/lib -lfann
else
    #stub file
    JBSNN=neo_jbs_nn
endif

OBJECTS = neo_globals.o \
          neo_nn_interface.o\
          neo_energy_grid.o \
          neo_interface.o \
          neo_profile_exp.o \
          neo_allocate_profile.o \
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
          neo_do.o \
          neo_error.o \
          neo_experimental_profiles.o \
          neo_init.o \
          neo_init_serial.o \
          neo_make_profiles.o \
          neo_map_experimental_profiles.o \
          neo_read_input.o \
          neo_harvest.o \
          neo_run.o \
          neo_spitzer.o \
          $(JBSNN).o

.SUFFIXES : .o .f90

all: neo_lib.a $(EXEC).o $(EXTRA_LIBS)
	$(FC) $(FFLAGS)  $(neo_libs) -o $(EXEC) $(EXEC).o neo_lib.a $(EXTRA_LIBS) $(LMATH) $(JBSNN_LIBS)

neo_lib.a: $(OBJECTS)
	$(ARCH) neo_lib.a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) $(neo_libs) -c $< $(JBSNN_LIBS)

clean:
	rm -rf *.o neo_lib.a $(EXEC) *~
	cd ${GACODE_ROOT}/modules ; rm -f neo*.mod

-include ../install/make.ext.${GACODE_PLATFORM}
