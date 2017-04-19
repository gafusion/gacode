include ${GACODE_ROOT}/shared/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC = gyro2gkv_exec

OBJECTS = gyro2gkv_globals.o \
          read_all.o \
          $(EXEC).o 

.SUFFIXES : .o .f90 

exec: $(OBJECTS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS) $(NETCDF)

.f90.o :
	$(FC) -I$(NETCDF_INC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(EXEC)
	rm -f ${GACODE_ROOT}/modules/gyro2gkv_globals.mod
