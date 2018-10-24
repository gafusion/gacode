include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

DESTLIB=${GACODE_ROOT}/python/gapy.so

SRC=	../src/gapy_realfluct.f90 \
	../src/gapy_torcut.f90 \
	../src/gapy_geo.f90 \
	../src/gapy_surfpar.f90 \
	../src/gapy_util.f90 \
	../src/gapy_vpro.f90

F2PYM:=$(filter-out --fcompiler=pg,$(F2PY))

ifndef F2PY
   all:
else
   all: objs
endif

objs:
	${F2PYM} --quiet -m gapy -h sgn.pyf ${SRC}
	${F2PY}  --quiet -lgomp --f90flags='${FOMP} ${FOPT} ${FMATH}' -c sgn.pyf ${SRC}

	rm -f ${DESTLIB}
	ln -s ${GACODE_ROOT}/gapy/f2py/gapy.so ${DESTLIB}

clean:
	rm -f gapy.so ${DESTLIB} sgn.pyf *~

exit:
	echo "No f2py"
