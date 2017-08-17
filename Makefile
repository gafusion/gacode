all:
	cd shared ; make
	cd le3 ; make
	cd gato ; make 
	cd neo ; make 
	cd vgen ; make 
	cd glf23 ; make
	cd tglf ; make
	cd cgyro ; make 
	cd gyro ; make 
	cd tgyro ; make
	cd profiles_gen ; make

clean:
	cd shared ; make clean
	cd le3 ; make clean
	cd gato ; make clean 
	cd neo ; make clean 
	cd vgen ; make clean 
	cd glf23 ; make clean
	cd tglf ; make clean
	cd cgyro ; make clean 
	cd gyro ; make clean 
	cd tgyro ; make clean
	cd profiles_gen ; make clean
	rm -f python/*/*.pyc
	rm -f python/*.pyc
	rm -f modules/*genmod*
	rm -f *.log
	rm -rf *regression_test/

distclean:
	cd shared ; make distclean
	cd neo ; make clean 
	cd tglf ; make clean
	cd cgyro ; make clean 
	cd gyro ; make distclean 
	cd tgyro ; make clean
	cd profiles_gen ; make clean
	rm -f python/*/*.pyc
	rm -f python/*.pyc
	rm -f modules/*genmod*
