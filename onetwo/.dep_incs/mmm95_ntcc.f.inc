mmm95_ntcc.o: mmm95_ntcc.f
	$(COMPILE) -I$(NETCDF_INC) -DONETWO mmm95_ntcc.f
