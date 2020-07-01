include ../platform/build/make.inc.${GACODE_PLATFORM}

all:
	cd math ; make "OPT=${OPT}"
	cd nclass ; make "OPT=${OPT}"
	cd UMFPACK ; make "OPT=${OPT}"
ifdef NEURAL_ROOT
	cd ${NEURAL_ROOT} && make "OPT=${OPT}"
endif

clean:
	cd math ; make clean
	cd nclass ; make clean
	cd UMFPACK ; make clean
ifdef NEURAL_ROOT
	cd ${NEURAL_ROOT} && make clean
endif

distclean:
	cd math ; make clean
