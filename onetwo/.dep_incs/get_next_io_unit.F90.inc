get_next_io_unit.o: error_handler.mod io_gcnmp.mod mpi_data.mod nrtype.mod get_next_io_unit.F90
	$(COMPILE) -I$(NETCDF_INC) -DONETWO get_next_io_unit.F90
