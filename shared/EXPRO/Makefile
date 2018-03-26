include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

LLIB=EXPRO_lib

OBJECTS = EXPRO_interface.o \
          EXPRO_locsim_interface.o \
          EXPRO_globals.o \
          EXPRO_alloc.o \
          EXPRO_alloc_control.o \
          EXPRO_compute_derived.o \
          EXPRO_locsim_alloc.o \
          EXPRO_locsim_profiles.o \
          EXPRO_palloc.o \
          EXPRO_pread.o \
          EXPRO_read.o \
          EXPRO_read_driver.o \
          EXPRO_read_geo.o \
          EXPRO_skip_header.o \
          EXPRO_write.o \
          EXPRO_write_derived.o \
          EXPRO_write_original.o

.SUFFIXES : .o .f90 

all: lib

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o *.lst $(LLIB).a ${GACODE_ROOT}/modules/expro_*.mod

