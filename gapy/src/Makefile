include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

DESTLIB=${GACODE_ROOT}/python/gapy.so

ifndef F2PY
   F2PY=f2py
endif

all:
	${F2PY} -m gapy -h sgn.pyf gapy_math.f90 gapy_fluct.f90
	${F2PY} -lgomp --f90flags='${FOMP}' -c sgn.pyf gapy_math.f90 gapy_fluct.f90 ${GACODE_ROOT}/shared/math/math_lib.a ${LMATH}
	rm -f ${DESTLIB}
	ln -s ${GACODE_ROOT}/gapy/src/gapy.so ${DESTLIB}

clean:
	rm -f gapy.so ${DESTLIB} sgn.pyf *~
