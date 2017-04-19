#----------------------------
# Toplevel makefile for cgyro
#----------------------------

include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = \
	${GACODE_ROOT}/shared/EXPRO/EXPRO_lib.a \
        ${GACODE_ROOT}/shared/GEO/GEO_lib.a \
        ${GACODE_ROOT}/shared/math/math_lib.a 

all:
	gacode_getversion > .VERSION
	cd ${GACODE_ROOT}/shared/math ; make
	cd ${GACODE_ROOT}/shared/GEO ; make
	cd ${GACODE_ROOT}/shared/EXPRO ; make
	cd src ; make

stub:
	cd src ; make stub

clean:
	cd ${GACODE_ROOT}/modules ; rm -f cgyro*.mod
	cd src ; make clean

.IGNORE:
