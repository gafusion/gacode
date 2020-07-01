include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG} -I$(NETCDF_INC)
else
   FFLAGS=${FOPT} -I$(NETCDF_INC)
endif

export EXTRA_LIBS = \
	${GACODE_ROOT}/le3/fluxfit/fluxfit_lib.a \
   	${GACODE_ROOT}/shared/math/math_lib.a 

EXEC = profiles_3d

OBJECTS = bound_deriv.o \
	cub_spline.o \
	profiles_3d.o

.SUFFIXES : .o .f90

all: $(EXEC)

$(EXEC): $(OBJECTS) $(EXTRA_LIBS)
	$(FC) $(FMATH) $(FFLAGS) -o $(EXEC) $(OBJECTS) $(EXTRA_LIBS) $(LMATH) $(NETCDF)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(EXEC)
