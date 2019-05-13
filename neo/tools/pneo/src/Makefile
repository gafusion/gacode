include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

#NEO_FLAGS = -cpp
ifdef CUDA_DIR
   NEO_FLAGS += -DNEO_HAVE_CUDA $(CUDA_INC)
   NEO_LIBS  += $(CUDA_LIB)
endif
ifdef PETSC_DIR
   NEO_FLAGS += -DNEO_HAVE_PETSC $(PETSC_INC)
   NEO_LIBS  += $(PETSC_LIB)
endif
ifdef SUPERLU_DIR
   NEO_FLAGS += -DNEO_HAVE_SUPERLU $(SUPERLU_INC)
   NEO_LIBS  += $(SUPERLU_LIB)
endif

EXEC1 = pneo
EXEC2 = pneo_rbf

OBJECTS = pneo_globals.o \
          $(EXEC1).o

.SUFFIXES : .o .f90

all: exec1 exec2

exec1: $(OBJECTS)
	$(FC) $(FFLAGS) $(NEO_FLAGS) -o $(EXEC1) $(OBJECTS) ${GACODE_ROOT}/neo/src/neo_lib.a $(EXTRA_LIBS) $(LMATH) $(NN_LIB) ${NEO_LIBS}


exec2: $(OBJECTS) $(EXEC2).o
	$(FC) $(FFLAGS) $(NEO_FLAGS) -o $(EXEC2) $(EXEC2).o $(LMATH)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) $(NEO_FLAGS) -c $<

clean:
	rm -f *.o $(EXEC1) $(EXEC2) ${GACODE_ROOT}/modules/pneo*.mod
