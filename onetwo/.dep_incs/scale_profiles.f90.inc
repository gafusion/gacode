scale_profiles.o: bc_values_gcnmp.mod common_constants.mod error_handler.mod grid_class.mod io_gcnmp.mod ions_gcnmp.mod mpi_data.mod nrtype.mod plasma_properties.mod tension_spline.mod vector_class.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/scale_profiles.f90
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/scale_profiles.f90
