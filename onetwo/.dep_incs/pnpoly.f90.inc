pnpoly.o: pnpoly.f90
	$(COMPILE) -I$(NETCDF_INC) -DONETWO pnpoly.f90
