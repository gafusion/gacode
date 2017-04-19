include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

LLIB=nclass_lib

OBJECTS = nclass_mod.o \
          rarray_copy.o \
          rarray_sum.o \
          rarray_zero.o \
          u_erf.o \
          u_lu_backsub.o \
          u_lu_decomp.o \
          write_mod.o 

.SUFFIXES : .o .f90 .f

all: lib

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f.o :
	$(F77) $(FMATH) $(FFLAGS) -c $<
 
clean:
	rm -rf *.o nclass_lib.a *~
