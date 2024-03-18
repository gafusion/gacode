USE_HARVEST=

include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

FC += ${FOMP} ${FACC}

EXEC = qlgyro

OBJECTS += qlgyro_globals.o \
	   qlgyro_allocate_globals.o \
	   qlgyro_cgyro_interface.o \
	   qlgyro_catch_error.o \
	   qlgyro_cgyro_cleanup.o \
	   qlgyro_cgyro_map.o \
	   qlgyro_comm_setup.o \
	   qlgyro_comm_sync.o \
	   qlgyro_ky_spectrum.o \
	   qlgyro_inout.o \
	   qlgyro_px0_spectrum.o \
	   qlgyro_read_input.o \
	   qlgyro_run_cgyro_balloon.o \
	   qlgyro_run_cgyro_fluxtube.o \
	   qlgyro_run.o \
	   qlgyro_sat1.o \
	   qlgyro_sat_mg.o \
	   qlgyro_sum_fluxes.o \
	   qlgyro_tglf_map.o \
	   qlgyro_trap_component_error.o

.SUFFIXES : .o .f90 .F90

all: qlgyro_lib.a $(EXEC).o $(EXTRA_LIBS)
	$(FC) $(FFLAGS) -o $(EXEC) $(EXEC).o qlgyro_lib.a $(EXTRA_LIBS) $(LMATH) $(NN_LIB) 

qlgyro_lib.a: $(OBJECTS)
	$(ARCH) qlgyro_lib.a $(OBJECTS)

.f90.o .F90.o:
	$(FC) $(FMATH) $(FFLAGS) $(NN_LIB) -c $<

clean:
	rm -f $(EXEC)
	rm -f *.o ${GACODE_ROOT}/modules/qlgyro*.mod
	rm -f qlgyro_lib.a
