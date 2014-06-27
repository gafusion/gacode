callglf2d.o: callglf2d.f
	$(COMPILE) -I$(NETCDF_INC) -DONETWO callglf2d.f
