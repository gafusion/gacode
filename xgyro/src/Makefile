include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

ifdef FFTW_INC
   FFLAGS += -I${FFTW_INC}
endif

FC += ${FOMP} ${FACC}

EXEC = xgyro

OBJECTS = xgyro_globals.o \
          xgyro_io.o \
          xgyro_mpi_setup.o \
          xgyro_read_input.o

.SUFFIXES : .o .f90 .F90

all: xgyro_lib.a $(EXEC).o $(EXTRA_LIBS)
	$(FC) $(FFLAGS) -o $(EXEC) $(EXEC).o xgyro_lib.a $(EXTRA_LIBS) $(LMATH)

xgyro_lib.a: $(OBJECTS)
	$(ARCH) xgyro_lib.a $(OBJECTS)

.f90.o .F90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -rf *.o xgyro_lib.a $(EXEC) *~ *.cub *.ptx *.gpu *.h *lst *.s

-include ../install/make.ext.${GACODE_PLATFORM}
