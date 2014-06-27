glf2d.o: glf2d.F
	$(COMPILE) -I$(NETCDF_INC) -DONETWO glf2d.F
