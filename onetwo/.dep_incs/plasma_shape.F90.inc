plasma_shape.o: common_constants.mod error_handler.mod io_gcnmp.mod io.mod mpi_data.mod nrtype.mod plasma_properties.mod vector_class.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/plasma_shape.F90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo plasma_shape.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/plasma_shape.F90
plasma_shape.mod: common_constants.mod error_handler.mod io_gcnmp.mod io.mod mpi_data.mod nrtype.mod plasma_properties.mod vector_class.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/plasma_shape.F90
