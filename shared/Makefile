include ../platform/build/make.inc.${GACODE_PLATFORM}

all:
	cd math ; make "OPT=${OPT}"
	cd GEO ; make "OPT=${OPT}"
	cd EXPRO ; make "OPT=${OPT}"
	cd fluxfit ; make "OPT=${OPT}"
	cd harvest_client ; make
ifdef FANN_ROOT
	cd neural; make
endif
	cd nclass ; make "OPT=${OPT}"
	cd UMFPACK ; make "OPT=${OPT}"

clean:
	cd math ; make clean
	cd GEO ; make clean
	cd EXPRO ; make clean
	cd fluxfit ; make clean
	cd harvest_client ; make clean
ifdef FANN_ROOT
	cd neural; make clean
endif
	cd nclass ; make clean
	cd UMFPACK ; make clean

distclean:
	cd math ; make clean
	cd GEO ; make clean
	cd EXPRO ; make clean
	cd harvest_client; make clean
