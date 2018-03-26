#----------------------------
# Toplevel makefile for vgen
#----------------------------

include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = \
	${GACODE_ROOT}/neo/src/neo_lib.a \
	${GACODE_ROOT}/shared/nclass/nclass_lib.a \
	${GACODE_ROOT}/shared/EXPRO/EXPRO_lib.a \
	${GACODE_ROOT}/shared/math/math_lib.a \
	${GACODE_ROOT}/shared/GEO/GEO_lib.a \
	${GACODE_ROOT}/shared/UMFPACK/UMFPACK_lib.a 

all:
	cd src ; make

clean:
	cd src ; make clean

.IGNORE:
