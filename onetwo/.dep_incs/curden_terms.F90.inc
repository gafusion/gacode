curden_terms.o: common_constants.mod dep_var.mod grid_class.mod ions_gcnmp.mod nrtype.mod numbrs.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/curden_terms.F90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo curden_terms.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/curden_terms.F90
curden_terms.mod: common_constants.mod dep_var.mod grid_class.mod ions_gcnmp.mod nrtype.mod numbrs.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/curden_terms.F90
