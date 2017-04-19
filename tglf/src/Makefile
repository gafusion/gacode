include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC1=tglf
EXEC2=tglf_mpi

LLIB=tglf_lib

ifdef FANN_ROOT
   # neural-net file
   TGLF_NN=tglf_nn_TM.fann
   TGLF_LIBS=-I$(GACODE_ROOT)/shared/neural/ -L$(GACODE_ROOT)/shared/neural/ -lbrainfuse -L${FANN_ROOT}/lib -lfann
else
   # stub file
   TGLF_NN=tglf_nn_TM
endif

OBJECTS = tglf_isnan.o      \
          tglf_modules.o    \
          tglf_pkg.o        \
          tglf_allocate.o   \
          tglf_deallocate.o \
          tglf_startup.o    \
          tglf_hermite.o    \
          tglf_inout.o      \
          tglf_setup_geometry.o \
          tglf_LS.o         \
          tglf_eigensolver.o \
          tglf_geometry.o   \
          tglf_matrix.o     \
          tglf_max.o        \
          tglf_interface.o  \
          tglf_error.o      \
          tglf_isinf.o      \
          tglf_shutdown.o   \
          tglf_read_input.o \
          tglf_multiscale_spectrum.o \
          tglf_kygrid.o     \
          tglf_run.o        \
          tglf_TM.o         \
          tglf_init_mpi.o   \
          tglf_run_mpi.o    \
          tglf_TM_mpi.o     \
          $(TGLF_NN).o

.SUFFIXES : .o .f90 .F90

all: $(LLIB).a $(EXEC1) $(EXEC2) Makefile

$(EXEC1): $(LLIB).a $(EXEC1).o
	$(FC) $(FFLAGS) -o $(EXEC1) $(EXEC1).o $(LLIB).a $(LMATH) $(TGLF_LIBS)
$(EXEC2): $(LLIB).a $(EXEC2).o
	$(FC) $(FFLAGS) -o $(EXEC2) $(EXEC2).o $(LLIB).a $(LMATH) $(TGLF_LIBS)

$(LLIB).a: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

tglf_run_mpi.o: tglf_run.F90
	$(FC) -o tglf_run_mpi.o $(FMATH) $(FFLAGS) -DMPI_TGLF=1 -c $<

.f90.o:
	$(FC) $(FMATH) $(FFLAGS) -c $< $(TGLF_LIBS)
.F90.o:
	$(FC) $(FMATH) $(FFLAGS) -c $< $(TGLF_LIBS)

clean:
	rm -f *.o *.lst $(LLIB).a $(EXEC1) $(EXEC2)
	rm -f ${GACODE_ROOT}/modules/tglf*.mod 

-include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}