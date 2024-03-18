#----------------------------
# Toplevel makefile for tgyro
#----------------------------

include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = \
       ${GACODE_ROOT}/neo/src/neo_lib.a \
       ${GACODE_ROOT}/shared/nclass/nclass_lib.a \
       ${GACODE_ROOT}/tglf/src/tglf_lib.a \
       ${GACODE_ROOT}/cgyro/src/cgyro_lib.a \
       ${GACODE_ROOT}/f2py/expro/expro_lib.a \
       ${GACODE_ROOT}/f2py/geo/geo_lib.a \
       ${GACODE_ROOT}/shared/math/math_lib.a \
       ${GACODE_ROOT}/shared/UMFPACK/UMFPACK_lib.a

all:
	gacode_getversion > .VERSION
	cd ${GACODE_ROOT}/shared ; make
	cd ${GACODE_ROOT}/f2py ; make
	cd ${GACODE_ROOT}/neo ; make
	cd ${GACODE_ROOT}/tglf ; make
	cd ${GACODE_ROOT}/cgyro ; make
	cd ${GACODE_ADD_ROOT}/qlgyro ; make
	cd src ; make

clean:
	cd src ; make clean
	rm -rf .VERSION

.IGNORE:
