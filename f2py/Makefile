#--------------------------------------
# Toplevel makefile for f2py components
#---------------------------------------

all:
	cd geo ; make
	cd expro ; make

pip:
	mkdir pygacode
	cp -r ${GACODE_ROOT}/python/cgyro pygacode
	cp -r ${GACODE_ROOT}/python/gyro pygacode
	cp -r ${GACODE_ROOT}/python/tgyro pygacode
	cp -r ${GACODE_ROOT}/python/neo pygacode
	cp -r ${GACODE_ROOT}/python/input_gacode pygacode
	cp -r ${GACODE_ROOT}/python/__init__.py pygacode

	python ./setup.py build_ext
	python ./setup.py sdist

so:
	cd module ; make

clean:
	cd geo ; make clean
	cd expro ; make clean
	cd module ; make clean
	cd test ; make clean
	rm -rf dist build MANIFEST sgn.pyf
	rm -rf pygacode

.IGNORE:
