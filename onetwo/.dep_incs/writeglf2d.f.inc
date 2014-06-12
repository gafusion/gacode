writeglf2d.o: writeglf2d.f
	$(COMPILE) -I$(NETCDF_INC) -DONETWO writeglf2d.f
