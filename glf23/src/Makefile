include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC=glf23
SRC=glf23

LLIB=glf23_lib

OBJECTS = glf23_module.o   \
          glf23_interface.o \
          glf23_init_mpi.o \
          glf23_isinf.o    \
          glf23_isnan_generic.o  \
          glf23_inout.o     \
          glf23_cgg.o       \
          glf23_TM.o         \
          glf23_run.o       \
          glf23_read_input.o 

.SUFFIXES : .o .f90 .F .f

all: $(LLIB).a $(EXEC) 

$(EXEC): $(LLIB).a $(SRC).o $(EXTRA_LIBS)
	$(FC) $(FFLAGS) -o $(EXEC) $(SRC).o $(LLIB).a $(EXTRA_LIBS) $(LMATH)

$(LLIB).a: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

.F.o :
	$(F77) $(FFLAGS)  -c $<

.f.o :
	$(F77) $(FFLAGS) -c $<

clean:
	rm -f *.o *.lst $(LLIB).a glf23
	rm -f ${GACODE_ROOT}/modules/glf23*.mod 

