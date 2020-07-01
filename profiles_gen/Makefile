include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

all:
	cd ${GACODE_ROOT}/f2py ; make
	cd locpargen ; make
	cd src ; make

clean:
	cd locpargen ; make clean
	cd src ; make clean
