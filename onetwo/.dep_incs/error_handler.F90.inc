error_handler.o: io_gcnmp.mod io.mod mpi_data.mod nrtype.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/error_handler.F90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo error_handler.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/error_handler.F90
error_handler.mod: io_gcnmp.mod io.mod mpi_data.mod nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/error_handler.F90
