source_terms_gcnmp.o: common_constants.mod grid_class.mod ions_gcnmp.mod nrtype.mod plasma_properties.mod solcon_gcnmp.mod vector_class.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/source_terms_gcnmp.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo source_terms_gcnmp.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/source_terms_gcnmp.f90
source_terms_gcnmp.mod: common_constants.mod grid_class.mod ions_gcnmp.mod nrtype.mod plasma_properties.mod solcon_gcnmp.mod vector_class.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/source_terms_gcnmp.f90
