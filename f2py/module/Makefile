include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

OBJ=gacode

SRC= 	../vis/vis.f90 \
	../geo/geo.f90 \
	../expro/expro.f90 \
	../expro/expro_util.f90 \
	../expro/expro_pycomm.f90

# Remove possible compilers
F2PYM:=$(filter-out --fcompiler=pg --fcompiler=intelem,$(F2PY))

ifndef F2PY
   all:
else
   all: objs
endif

objs:
	${F2PYM} --quiet --overwrite-signature -m ${OBJ} -h sgn.pyf ${SRC}
	${F2PY}  --quiet -lgomp -c sgn.pyf ${SRC}

	cp *.so ${GACODE_ROOT}/f2py/pygacode/

clean:
	rm -f sgn.pyf *~ *.so ${GACODE_ROOT}/f2py/pygacode/*.so

exit:
	echo "No f2py"
