include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

DESTLIB=${GACODE_ROOT}/python/gapy.so

LMATH:=$(filter-out -lfftw3,$(LMATH))
LMATH:=$(filter-out -lfftw3f,$(LMATH))
LMATH:=$(filter-out -lfftw3_threads,$(LMATH))
LMATH:=$(filter-out -lfftw3f_threads,$(LMATH))
LMATH:=$(filter-out -lcufft,$(LMATH))
F2PYM:=$(filter-out --fcompiler=pg,$(F2PY))

ifndef F2PY
   all:
else
   all: objs
endif

objs:
	${F2PYM} --quiet -m gapy -h sgn.pyf gapy_math.f90 gapy_fluct.f90
	${F2PY} --quiet -lgomp --f90flags='${FOMP} ${FOPT}' -c sgn.pyf gapy_math.f90 gapy_fluct.f90 ${GACODE_ROOT}/shared/math/math_lib.a ${LMATH}

	rm -f ${DESTLIB}
	ln -s ${GACODE_ROOT}/gapy/src/gapy.so ${DESTLIB}

clean:
	rm -f gapy.so ${DESTLIB} sgn.pyf *~

exit:
	echo "No f2py"
