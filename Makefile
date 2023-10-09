all:
	printf "GACODE BRANCH: %s\n" `git branch | cut -d' ' -f2` > .VERSION
	gacode_getversion >> .VERSION
#	printf "GACODE HASH: %s\n" `git rev-parse \`git branch | cut -d' ' -f2\`` >> GACODE_BUILD_INFO
#	printf 'COMPILE DATE: %(%a %b %e %H:%M:%S %Z %Y)T\n' >> GACODE_BUILD_INFO
	cd shared ; make
#	cd le3 ; make
	cd neo ; make
	cd vgen ; make
	cd glf23 ; make
	cd tglf ; make
	cd cgyro ; make
	cd gyro ; make
	cd tgyro ; make
	cd profiles_gen ; make
	cd f2py ; make
	@echo "GACODE build done"

clean:
	rm -f .VERSION
	cd shared ; make clean
#	cd le3 ; make clean
	cd neo ; make clean
	cd vgen ; make clean
	cd glf23 ; make clean
	cd tglf ; make clean
	cd cgyro ; make clean
	cd gyro ; make clean
	cd tgyro ; make clean
	cd profiles_gen ; make clean
	cd f2py ; make clean
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

distclean:
	cd shared ; make distclean
	cd neo ; make clean
	cd tglf ; make clean
	cd cgyro ; make clean
	cd gyro ; make clean
	cd tgyro ; make clean
	cd profiles_gen ; make clean
	cd f2py ; make clean
	rm -f python/*/*.pyc
	rm -f python/*.pyc
	rm -f modules/*genmod*
	rm -f modules/*.mod
