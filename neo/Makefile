#---------------------------
# Toplevel makefile for neo
#---------------------------

include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = \
       ${GACODE_ROOT}/shared/EXPRO/EXPRO_lib.a \
       ${GACODE_ROOT}/shared/math/math_lib.a \
       ${GACODE_ROOT}/shared/GEO/GEO_lib.a \
       ${GACODE_ROOT}/shared/UMFPACK/UMFPACK_lib.a \
       ${GACODE_ROOT}/shared/nclass/nclass_lib.a

all:
	gacode_getversion > .VERSION
	cd ${GACODE_ROOT}/shared ; make
	cd src ; make
	cd ${GACODE_ROOT}/neo/tools/pneo ; make

clean:
	cd src ; make clean
	cd ${GACODE_ROOT}/neo/tools/pneo ; make clean

.IGNORE:
