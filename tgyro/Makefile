#----------------------------
# Toplevel makefile for tgyro
#----------------------------

include ../platform/build/make.inc.${GACODE_PLATFORM}
include Makefile.extra

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
