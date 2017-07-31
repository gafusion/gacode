#----------------------------
# Toplevel makefile for vgen
#----------------------------

include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = \
	${GACODE_ROOT}/cgyro/src/cgyro_lib.a \
	${GACODE_ROOT}/shared/GEO/GEO_lib.a \
	${GACODE_ROOT}/shared/EXPRO/EXPRO_lib.a \
	${GACODE_ROOT}/shared/math/math_lib.a

all:
	cd src ; make

clean:
	cd src ; make clean

.IGNORE:
