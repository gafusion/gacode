all:
	gacode_getversion > .VERSION
	cd shared ; make
	cd neo ; make
	cd vgen ; make
	cd tglf ; make
	cd cgyro ; make
	cd xgyro ; make
	cd tgyro ; make
	cd profiles_gen ; make
	cd f2py ; make
	cd qlgyro ; make

# Clean but do not build gyro and le3

clean:
	rm -f .VERSION
	cd shared ; make clean
	cd shared/hybridtest ; make clean
	cd neo ; make clean
	cd vgen ; make clean
	cd tglf ; make clean
	cd cgyro ; make clean
	cd xgyro ; make clean
	cd tgyro ; make clean
	cd profiles_gen ; make clean
	cd f2py ; make clean
	cd qlgyro ; make clean
	cd gyro ; make clean
	cd le3 ; make clean
	rm -f f2py/*/*.pyc
	rm -f f2py/*.pyc
	rm -f f2py/pygacode/*/*.pyc
	rm -rf f2py/*/__pycache__
	rm -rf f2py/__pycache__
	rm -rf f2py/pygacode/*/__pycache__
	rm -f modules/*genmod*
	rm -f *.log
	rm -rf *regression_test/
	rm -rf python
