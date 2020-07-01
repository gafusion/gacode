include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC = ptglf

OBJECTS = ptglf_globals.o \
          ptglf_init.o \
          $(EXEC).o

.SUFFIXES : .o .f90

all: exec

exec: $(OBJECTS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS) ${GACODE_ROOT}/tglf/src/tglf_lib.a $(EXTRA_LIBS) $(LMATH)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(EXEC) ${GACODE_ROOT}/modules/ptglf*.mod *~
