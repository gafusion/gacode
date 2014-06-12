fusion_gcnmp.o: common_constants.mod dep_var.mod fast_ion_data_gcnmp.mod grid_class.mod ions_gcnmp.mod nrtype.mod solcon_gcnmp.mod source_terms_gcnmp.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/fusion_gcnmp.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo fusion_gcnmp.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/fusion_gcnmp.f90
fusion_gcnmp.mod: common_constants.mod dep_var.mod fast_ion_data_gcnmp.mod grid_class.mod ions_gcnmp.mod nrtype.mod solcon_gcnmp.mod source_terms_gcnmp.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/fusion_gcnmp.f90
