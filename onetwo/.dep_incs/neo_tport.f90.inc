neo_tport.o: common_constants.mod curden_terms.mod dep_var.mod error_handler.mod grid_class.mod io_gcnmp.mod ions_gcnmp.mod mpi_data.mod neutral_beams.mod nrtype.mod plasma_properties.mod solcon_gcnmp.mod source_terms_gcnmp.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/neo_tport.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo neo_tport.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/neo_tport.f90
neo_tport.mod: common_constants.mod curden_terms.mod dep_var.mod error_handler.mod grid_class.mod io_gcnmp.mod ions_gcnmp.mod mpi_data.mod neutral_beams.mod nrtype.mod plasma_properties.mod solcon_gcnmp.mod source_terms_gcnmp.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/neo_tport.f90
