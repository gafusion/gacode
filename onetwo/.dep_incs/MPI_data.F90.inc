MPI_data.o: common_constants.mod nrtype.mod solcon_gcnmp.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/MPI_data.F90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo mpi_data.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/MPI_data.F90
mpi_data.mod: common_constants.mod nrtype.mod solcon_gcnmp.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/MPI_data.F90
