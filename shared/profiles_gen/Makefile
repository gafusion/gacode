include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

all:
	cd ${GACODE_ROOT}/shared/math ; make
	cd ${GACODE_ROOT}/shared/GEO ; make
	cd ${GACODE_ROOT}/shared/EXPRO ; make
	cd ${GACODE_ROOT}/shared/fluxfit ; make
	cd locpargen ; make
	cd src ; make
	cd exprotools ; make

clean:
	cd locpargen ; make clean
	cd src ; make clean
	cd exprotools ; make clean
