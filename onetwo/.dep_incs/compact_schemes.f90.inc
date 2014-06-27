compact_schemes.o: common_constants.mod error_handler.mod io_gcnmp.mod mpi_data.mod nrtype.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/compact_schemes.f90
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/compact_schemes.f90
