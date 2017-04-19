include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

LLIB=GEO_lib

OBJECTS = GEO_interface.o \
          GEO_alloc.o \
          GEO_do.o \
          GEO_interp.o \
          GEO_write.o

.SUFFIXES : .o .f90 

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(LLIB).a ${GACODE_ROOT}/modules/geo_interface.mod
