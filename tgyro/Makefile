#----------------------------
# Toplevel makefile for tgyro
#----------------------------

include ../platform/build/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = \
       ${GACODE_ROOT}/neo/src/neo_lib.a \
       ${GACODE_ROOT}/shared/nclass/nclass_lib.a \
       ${GACODE_ROOT}/tglf/src/tglf_lib.a \
       ${GACODE_ROOT}/glf23/src/glf23_lib.a \
       ${GACODE_ROOT}/gyro/src/gyro_lib.a \
       ${GACODE_ROOT}/gyro/BLEND/BLEND_lib.a \
       ${GACODE_ROOT}/cgyro/src/cgyro_lib.a \
       ${GACODE_ROOT}/gyro/TRANSPOSE/TRANSP_lib.a \
       ${GACODE_ROOT}/gyro/GKEIGEN/GKEIGEN_lib.a \
       ${GACODE_ROOT}/gyro/SSUB/SSUB_lib.a \
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
	cd ${GACODE_ROOT}/glf23 ; make
	cd ${GACODE_ROOT}/gyro ; make
	cd ${GACODE_ROOT}/cgyro ; make
	cd src ; make

clean:
	cd src ; make clean
	rm -rf .VERSION

.IGNORE:
