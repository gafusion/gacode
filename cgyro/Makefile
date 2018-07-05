#----------------------------
# Toplevel makefile for cgyro
#----------------------------

include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = \
	${GACODE_ROOT}/gapy/EXPRO/EXPRO_lib.a \
        ${GACODE_ROOT}/gapy/src/geo_lib.a \
        ${GACODE_ROOT}/shared/math/math_lib.a

all:
	gacode_getversion > .VERSION
	cd ${GACODE_ROOT}/shared/math ; make
	cd ${GACODE_ROOT}/gapy ; make
	cd src ; make
	cd ${GACODE_ROOT}/cgyro/tools/globalmap ; make

clean:
	cd ${GACODE_ROOT}/modules ; rm -f cgyro*.mod
	cd src ; make clean
	cd ${GACODE_ROOT}/cgyro/tools/globalmap ; make clean

deepclean:
	cd ${GACODE_ROOT}/shared/math ; make clean
	cd ${GACODE_ROOT}/gapy ; make clean
	cd ${GACODE_ROOT}/modules ; rm -f cgyro*.mod
	cd src ; make clean
	cd ${GACODE_ROOT}/cgyro/tools/globalmap ; make clean

.IGNORE:
