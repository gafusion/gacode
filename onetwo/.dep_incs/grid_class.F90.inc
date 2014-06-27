grid_class.o: common_constants.mod error_handler.mod io_gcnmp.mod io.mod mpi_data.mod nrtype.mod plasma_properties.mod solcon_gcnmp.mod tension_spline.mod vector_class.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/grid_class.F90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo grid_class.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/grid_class.F90
grid_class.mod: common_constants.mod error_handler.mod io_gcnmp.mod io.mod mpi_data.mod nrtype.mod plasma_properties.mod solcon_gcnmp.mod tension_spline.mod vector_class.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/grid_class.F90
