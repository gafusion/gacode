#--------------------------------------
# Toplevel makefile for f2py components
#---------------------------------------

all:
	cd geo ; make
	cd expro ; make

vismod:
	cd vis ; make

pip:
	python setup.py build_ext

clean:
	cd geo ; make clean
	cd expro ; make clean
	cd vis ; make clean
	rm -rf dist build *egg-info *.so

.IGNORE:
