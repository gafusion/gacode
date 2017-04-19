include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

LLIB=BLEND_lib

OBJECTS = BLEND_private.o \
          BLEND_init.o \
          BLEND_do.o \
          BLEND_F.o \
          BLEND_Fp.o \
          BLEND_f2.o \
          BLEND_f2p.o \
          BLEND_f3.o \
          BLEND_f3p.o \
          BLEND_f4.o \
          BLEND_f4p.o \
          BLEND_cleanup.o 

.SUFFIXES : .o .f90 

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -rf *.o $(LLIB).a 
	rm -rf ${GACODE_ROOT}/modules/blend*.mod
