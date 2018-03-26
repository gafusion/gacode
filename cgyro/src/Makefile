include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

ifdef FFTW_INC
   FFLAGS += -I${FFTW_INC}
endif

FC += ${FOMP}

EXEC = cgyro

OBJECTS = cgyro_globals.o \
          cgyro_parallel_lib.o \
          cgyro_io.o \
          cgyro_advect_wavenumber.o \
	  cgyro_init_rotation.o \
          cgyro_equilibrium.o \
          cgyro_check.o \
	  cgyro_check_memory.o \
          cgyro_error_estimate.o \
          cgyro_field.o \
          cgyro_field_coefficients.o \
          cgyro_flux.o \
          cgyro_freq.o \
          cgyro_upwind.o \
          cgyro_init.o \
          cgyro_init_h.o \
          cgyro_init_arrays.o \
	  cgyro_init_collision_simple.o \
          cgyro_init_collision.o \
          cgyro_init_manager.o \
	  cgyro_step_collision_simple.o \
          cgyro_step_collision.o \
          cgyro_kernel.o \
          cgyro_mpi_grid.o \
	  cgyro_make_profiles.o \
	  cgyro_nl_direct.o \
	  cgyro_nl_fftw.o \
          cgyro_read_input.o \
          cgyro_rhs.o \
          cgyro_run.o \
          cgyro_shear.o \
          cgyro_shear_interface.o \
          cgyro_shear_hammett.o \
	  cgyro_step_gk.o \
          cgyro_write_initdata.o \
          cgyro_read_restart.o \
          cgyro_write_restart.o \
          cgyro_write_timedata.o \
          cgyro_write_hosts.o

.SUFFIXES : .o .f90 .F90

all: cgyro_lib.a $(EXEC).o $(EXTRA_LIBS)
	$(FC) $(FFLAGS) -o $(EXEC) $(EXEC).o cgyro_lib.a $(EXTRA_LIBS) $(LMATH)

cgyro_lib.a: $(OBJECTS)
	$(ARCH) cgyro_lib.a $(OBJECTS)

.f90.o .F90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -rf *.o cgyro_lib.a $(EXEC) *~ *.cub *.ptx *.gpu *.h *lst

-include ../install/make.ext.${GACODE_PLATFORM}
