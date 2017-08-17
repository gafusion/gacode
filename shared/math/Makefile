include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

LLIB=math_lib

OBJECTS = math_constants.o \
          bound_deriv.o \
          bound_extrap.o \
          cub_spline.o \
          cub_spline_deriv.o \
          cub_spline_deriv2.o \
          energy_integral.o \
          gacode_system.o \
          gauss_legendre.o \
          gauss_integ.o \
          i0.o \
          math_scaleint.o \
          math_scaleintv.o \
          parallel_dim.o \
          polydiff.o \
          poly2diff.o \
          proc_time.o \
          pseudospec_lib.o \
          timer_lib.o

.SUFFIXES : .o .f90 .f

all: lib

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f.o :
	$(F77) $(FFLAGS) -c $<

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(LLIB).a *.mod
	rm -f ${GACODE_ROOT}/modules/math_constants.mod
	rm -f ${GACODE_ROOT}/modules/timer_lib.mod

-include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}
