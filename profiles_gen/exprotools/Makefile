include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC1 = exprodump
EXEC2 = create

EXTRA_LIBS = ${GACODE_ROOT}/shared/EXPRO/EXPRO_lib.a \
             ${GACODE_ROOT}/shared/GEO/GEO_lib.a  \
             ${GACODE_ROOT}/shared/math/math_lib.a 

OBJECTS = create_globals.o \
          create_check.o \
          create_importvec.o \
          create_set.o \
          create_read_input.o \
          create_write.o 

.SUFFIXES : .o .f90

all : $(EXEC1) $(EXEC2)

$(EXEC1): $(EXEC1).o $(OBJECTS) $(EXTRA_LIBS)
	$(FC) $(FFLAGS) -o $(EXEC1) $(EXEC1).o $(EXTRA_LIBS) $(LMATH)

$(EXEC2): $(EXEC2).o $(OBJECTS) $(EXTRA_LIBS)
	$(FC) $(FFLAGS) -o $(EXEC2) $(EXEC2).o $(OBJECTS) $(EXTRA_LIBS) $(LMATH)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(EXEC1) $(EXEC2) ${GACODE_ROOT}/modules/create*.mod
