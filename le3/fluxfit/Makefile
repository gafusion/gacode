include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

LLIB=fluxfit_lib

OBJECTS = fluxfit_globals.o \
          fluxfit_driver.o \
          fluxfit_error.o \
          fluxfit_f_model.o \
          fluxfit_fourier.o \
          fluxfit_intersect.o \
          fluxfit_minmax.o \
          fluxfit_moments.o 

.SUFFIXES : .o .f90 

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(LLIB).a ${GACODE_ROOT}/modules/fluxfit*.mod
