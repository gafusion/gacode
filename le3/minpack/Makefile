include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

LLIB=minpack_lib

OBJECTS = dogleg.o \
	dpmpar.o \
        enorm.o \
        fdjac1.o \
        hybrd1.o \
        hybrd.o \
        qform.o \
        qrfac.o \
        r1mpyq.o \
        r1updt.o 

.SUFFIXES : .o .f

all: lib

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f.o :
	$(F77) $(FFLAGS) -c $<
 
clean:
	rm -rf *.o $(LLIB).a *~
