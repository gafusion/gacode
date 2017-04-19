include ../platform/build/make.inc.${GACODE_PLATFORM}

-include install/make.ext.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG} ${GATOFLAGS}
else
   FFLAGS=${FOPT} ${GATOFLAGS}
endif

EXEC = gato

OBJECTS = $(EXEC).o gato_plotg_library.stub.o

.SUFFIXES : .o .f

all: $(OBJECTS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS)

.f.o :
	$(F77) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o *~ $(EXEC)

