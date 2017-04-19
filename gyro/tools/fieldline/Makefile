include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC = poincare

OBJ = fieldline_input_data.o \
      fieldline_runge_kutta.o \
      func.o \
      read_input.o \
      setup.o

.SUFFIXES : .o .f90 

all: $(OBJ) $(EXEC).o
	$(FC) $(FFLAGS) -o $(EXEC) $(EXEC).o $(OBJ)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:  
	rm -f *.o $(EXEC)
	rm -f ${GACODE_ROOT}/modules/fieldline*.mod
