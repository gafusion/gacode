get_nboutput.o: get_nboutput.f90
	$(COMPILE) -I$(NETCDF_INC) -DONETWO get_nboutput.f90
