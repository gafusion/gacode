include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG} -I$(NETCDF_INC)
else
   FFLAGS=${FOPT} -I$(NETCDF_INC)
endif

export EXTRA_LIBS = \
	${GACODE_ROOT}/shared/math/math_lib.a \
	${GACODE_ROOT}/shared/EXPRO/EXPRO_lib.a \
	${GACODE_ROOT}/shared/GEO/GEO_lib.a \
	${GACODE_ROOT}/shared/fluxfit/fluxfit_lib.a

EXEC = prgen

OBJECTS = prgen_globals.o \
          allocate_internals.o \
          allocate_iterdb_vars.o \
          allocate_plasmastate_vars.o \
	  allocate_peqdsk_vars.o \
	  allocate_corsica_vars.o \
	  allocate_ufile_vars.o \
          prgen_get_chi.o \
          prgen_ion_name.o \
          prgen_map_iterdb.o \
          prgen_map_null.o \
          prgen_map_plasmastate.o \
	  prgen_map_peqdsk.o \
	  prgen_map_corsica.o \
	  prgen_map_ufile.o \
	  prgen_map_inputprofiles.o \
          prgen_read_cer.o \
          prgen_read_dskgato.o \
          prgen_read_omfit.o \
          prgen_read_gato.o \
          prgen_read_iterdb.o \
          prgen_read_iterdb_nc.o \
          prgen_read_null.o \
          prgen_read_plasmastate.o \
	  prgen_read_peqdsk.o \
	  prgen_read_corsica.o \
	  prgen_read_ufile.o \
	  prgen_read_inputprofiles.o \
          prgen_write.o 

.SUFFIXES : .o .f90

all: $(EXEC)

$(EXEC): $(OBJECTS) $(EXEC).o $(EXTRA_LIBS)
	$(FC) $(FMATH) $(FFLAGS) -o $(EXEC) $(EXEC).o $(OBJECTS) $(EXTRA_LIBS) $(LMATH) $(NETCDF)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(EXEC)
	rm -f ${GACODE_ROOT}/modules/prgen*.mod
