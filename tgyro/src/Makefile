include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

FC += ${FOMP}

EXEC = tgyro_main

ifdef FANN_ROOT
   TGYRO_NN=tgyro_eped_nn.fann
   TGYRO_LIBS=-L$(GACODE_ROOT)/shared/harvest_client/ \
              -I$(GACODE_ROOT)/shared/harvest_client/ \
              -L$(GACODE_ROOT)/shared/neural/ \
              -I$(GACODE_ROOT)/shared/neural/ \
              -lbrainfuse -lharvest \
              -L$(FANN_ROOT)/lib -lfann
else
   TGYRO_NN=tgyro_eped_nn
   TGYRO_LIBS=-L$(GACODE_ROOT)/shared/harvest_client/ \
              -I$(GACODE_ROOT)/shared/harvest_client/ \
              -lharvest
endif

OBJECTS += tgyro_globals.o \
           tgyro_iteration_variables.o \
           tgyro_allocate_globals.o \
	   tgyro_ped.o \
           tgyro_auxiliary_routines.o \
           tgyro_catch_error.o \
           tgyro_comm_setup.o \
           tgyro_comm_sync.o \
           tgyro_expro_map.o \
           tgyro_flux.o \
           tgyro_flux_set.o \
           tgyro_flux_vector.o \
           tgyro_gyro_map.o \
           tgyro_harvest.o \
           tgyro_ifs_pppl.o \
           tgyro_init_profiles.o \
           tgyro_init_profiles_axis.o \
           tgyro_iteration_driver.o \
           tgyro_iteration_standard.o \
           tgyro_iteration_serial.o \
           tgyro_iteration_parallel.o \
           tgyro_iteration_simplerelax.o \
           tgyro_multi_driver.o \
           tgyro_neo_map.o \
           tgyro_profile_functions.o \
           tgyro_profile_reintegrate.o \
           tgyro_profile_set.o \
           tgyro_pressure.o \
           tgyro_quasigrad.o \
           tgyro_read_input.o \
           tgyro_read_routines.o \
           tgyro_residual.o \
           tgyro_restart.o \
           tgyro_source.o \
           tgyro_target_vector.o \
           tgyro_tglf_map.o \
           tgyro_glf23_map.o \
           tgyro_trap_component_error.o \
           tgyro_volume_int.o \
           tgyro_write_input.o \
           tgyro_write_intermediate.o \
           tgyro_write_data.o \
           $(TGYRO_NN).o \
           $(EXEC).o

.SUFFIXES : .o .f90

all: $(OBJECTS) $(EXTRA_LIBS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS) $(EXTRA_LIBS) $(LMATH) $(TGYRO_LIBS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) $(TGYRO_LIBS) -c $<  

clean:
	rm -f $(EXEC)
	rm -f *.o ${GACODE_ROOT}/modules/tgyro*.mod

