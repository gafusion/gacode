dep_var_gcnmp.o: common_constants.mod nrtype.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/dep_var_gcnmp.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo dep_var.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/dep_var_gcnmp.f90
dep_var.mod: common_constants.mod nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/dep_var_gcnmp.f90
