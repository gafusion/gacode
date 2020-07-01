#----------------------------
# Toplevel makefile for cgyro
#----------------------------

include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = \
	${GACODE_ROOT}/f2py/expro/expro_lib.a \
        ${GACODE_ROOT}/f2py/geo/geo_lib.a \
        ${GACODE_ROOT}/shared/math/math_lib.a

all:
	gacode_getversion > .VERSION
	cd ${GACODE_ROOT}/shared/math ; make
	cd ${GACODE_ROOT}/f2py ; make
	cd src ; make

clean:
	cd ${GACODE_ROOT}/modules ; rm -f cgyro*.mod
	cd src ; make clean

deepclean:
	cd ${GACODE_ROOT}/shared/math ; make clean
	cd ${GACODE_ROOT}/f2py ; make clean
	cd ${GACODE_ROOT}/modules ; rm -f cgyro*.mod
	cd src ; make clean

.IGNORE:
