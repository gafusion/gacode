include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}
include ../Makefile.extra

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

FC += ${FOMP} ${FACC}

EXEC = tgyro_main

ifdef FANN_ROOT
   NN=tgyro_eped_nn.fann
else
   NN=tgyro_eped_nn
endif

ifdef MMM_HOME
  MMM_OBJ += tgyro_mmm_map.o
  MMM_LIB = -L$(MMM_HOME)/lib -lmmm
  FFLAGS += -DMMM -I$(MMM_HOME)/include
else
  OBJECTS += tgyro_mmm_dummy.o
endif

ifdef GACODE_ADD_ROOT
  EXTRA_LIBS += $(GACODE_ADD_ROOT)/qlgyro/src/qlgyro_lib.a
endif

OBJECTS += tgyro_globals.o \
           tgyro_iteration_variables.o \
           tgyro_allocate_globals.o \
	   tgyro_ped.o \
           tgyro_auxiliary_routines.o \
           tgyro_catch_error.o \
	   tgyro_cgyro_map.o \
           tgyro_comm_setup.o \
           tgyro_comm_sync.o \
           tgyro_etgcrit.o \
           tgyro_expro_map.o \
           $(MMM_OBJ) \
           tgyro_flux.o \
           tgyro_flux_set.o \
           tgyro_flux_vector.o \
           tgyro_ifs_pppl.o \
           tgyro_init_profiles.o \
           tgyro_init_profiles_axis.o \
           tgyro_iteration_driver.o \
           tgyro_iteration_standard.o \
           tgyro_iteration_diagonal.o \
           tgyro_iteration_serial.o \
           tgyro_iteration_parallel.o \
           tgyro_iteration_simplerelax.o \
           tgyro_iteration_zero.o \
           tgyro_multi_driver.o \
	   tgyro_cgyro_iterate.o \
           tgyro_neo_map.o \
           tgyro_profile_functions.o \
           tgyro_profile_reintegrate.o \
           tgyro_profile_set.o \
           tgyro_pressure.o \
           tgyro_quasigrad.o \
           tgyro_rad.o \
           tgyro_read_input.o \
           tgyro_read_routines.o \
           tgyro_residual.o \
           tgyro_restart.o \
           tgyro_source.o \
           tgyro_target_vector.o \
           tgyro_tglf_map.o \
           tgyro_trap_component_error.o \
           tgyro_volume_int.o \
           tgyro_write_input.o \
           tgyro_write_intermediate.o \
           tgyro_write_data.o \
           $(NN).o \
           $(EXEC).o

.SUFFIXES : .o .f90

all: $(OBJECTS) $(EXTRA_LIBS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS) $(EXTRA_LIBS) $(LMATH) $(NN_LIB) $(MMM_LIB)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) $(NN_LIB) -c $<

clean:
	rm -f $(EXEC)
	rm -f *.o ${GACODE_ROOT}/modules/tgyro*.mod

