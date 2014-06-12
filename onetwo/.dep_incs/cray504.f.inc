cray504.o: param.mod cray504.f
	$(COMPILE) -I$(NETCDF_INC) -DONETWO cray504.f
