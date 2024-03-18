include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG} -I$(NETCDF_INC)
else
   FFLAGS=${FOPT} -I$(NETCDF_INC)
endif

export EXTRA_LIBS = \
	${GACODE_ROOT}/f2py/expro/expro_lib.a \
        ${GACODE_ROOT}/f2py/geo/geo_lib.a

EXEC = prgen

OBJECTS = prgen_globals.o \
          prgen_allocate.o \
          prgen_get_chi.o \
          prgen_geometry.o \
          prgen_ion_name.o \
          prgen_map_iterdb.o \
          prgen_map_manual.o \
          prgen_map_null.o \
          prgen_map_plasmastate.o \
	  prgen_map_peqdsk.o \
	  prgen_map_genf.o \
	  prgen_map_inputgacode.o \
          prgen_read_cer.o \
          prgen_read_ti.o \
          prgen_read_iterdb.o \
          prgen_read_iterdb_nc.o \
          prgen_read_manual.o \
          prgen_read_plasmastate.o \
	  prgen_read_peqdsk.o \
	  prgen_read_genf.o \
	  prgen_read_inputgacode.o \
	  prgen_read_inputprofiles.o \
          prgen_swap.o \
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
