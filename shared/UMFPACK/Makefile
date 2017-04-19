include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

# We don't apply debug flags to UMFPACK because its old and trusted.

ifeq ($(OPT),debug)
   FFLAGS=
else
   FFLAGS=${FOPT}
endif
      
LLIB = UMFPACK_lib

HARWELL	= mc21b.o mc13e.o

UMFD = umd2fb.o \
       umd2co.o \
       umd2fa.o \
       umd2f0.o \
       umd2f1.o \
       umd2f2.o \
       umd2fg.o \
       umd2in.o \
       umd2of.o \
       umd2s2.o \
       umd2sl.o \
       umd2so.o \
       umd2su.o \
       umd2er.o \
       umd2p1.o \
       umd2p2.o \
       umd2lt.o \
       umd2ut.o \
       umd2rf.o \
       umd2ra.o \
       umd2r0.o \
       umd2r2.o \
       umd2rg.o \
       umd21i.o \
       umz2fb.o \
       umz2co.o \
       umz2fa.o \
       umz2f0.o \
       umz2f1.o \
       umz2f2.o \
       umz2fg.o \
       umz2in.o \
       umz2of.o \
       umz2s2.o \
       umz2sl.o \
       umz2so.o \
       umz2su.o \
       umz2er.o \
       umz2p1.o \
       umz2p2.o \
       umz2rf.o \
       umz2ra.o \
       umz2r0.o \
       umz2r2.o \
       umz2rg.o \
       umz21i.o

.SUFFIXES: .o .f

lib: $(UMFD) $(HARWELL)
	$(ARCH) $(LLIB).a $(UMFD) $(HARWELL)

.f.o:
	$(F77) $(FFLAGS) -c $<

clean:
	rm -f *.o $(LLIB).a
