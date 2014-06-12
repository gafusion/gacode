get_tot_pressure.o: common_constants.mod dep_var.mod error_handler.mod fast_ion_data_gcnmp.mod grid_class.mod io_gcnmp.mod ions_gcnmp.mod neutral_beams.mod nrtype.mod plasma_properties.mod vector_class.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/get_tot_pressure.f90
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/get_tot_pressure.f90
