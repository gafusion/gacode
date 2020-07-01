include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = \
       ${GACODE_ROOT}/le3/minpack/minpack_lib.a \
       ${GACODE_ROOT}/f2py/geo/geo_lib.a \
       ${GACODE_ROOT}/le3/fluxfit/fluxfit_lib.a

all:
	gacode_getversion > .VERSION
	cd ${GACODE_ROOT}/shared/math ; make
	cd ${GACODE_ROOT}/f2py ; make
	cd ${GACODE_ROOT}/le3/fluxfit ; make
	cd minpack ; make
	cd src ; make
	cd profiles_3d ; make

clean:
	cd minpack ; make clean
	cd src ; make clean
	cd profiles_3d ; make clean
	cd ../modules ; rm -f le3*.mod
	rm -f .VERSION
.IGNORE:
