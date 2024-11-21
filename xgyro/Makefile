#----------------------------
# Toplevel makefile for xgyro
#----------------------------

include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = \
	${GACODE_ROOT}/cgyro/src/cgyro_lib.a \
	${GACODE_ROOT}/f2py/expro/expro_lib.a \
        ${GACODE_ROOT}/f2py/geo/geo_lib.a \
        ${GACODE_ROOT}/shared/math/math_lib.a

all:
	gacode_getversion > .VERSION
	cd ${GACODE_ROOT}/cgyro && make
	cd src && make

clean:
	cd ${GACODE_ROOT}/modules && rm -f xgyro*.mod
	cd src && make clean

deepclean:
	cd ${GACODE_ROOT}/cgyro && make deepclean
	cd ${GACODE_ROOT}/modules && rm -f xgyro*.mod
	cd src && make clean

.IGNORE:
