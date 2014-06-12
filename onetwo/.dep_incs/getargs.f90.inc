getargs.o: error_handler.mod gcnmp_version.mod io_gcnmp.mod nrtype.mod getargs.f90
	$(COMPILE) -I$(NETCDF_INC) -DONETWO getargs.f90
