thermalization.o: common_constants.mod constnts.mod error_handler.mod io_gcnmp.mod nrtype.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/thermalization.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo thermalization.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/thermalization.f90
thermalization.mod: common_constants.mod constnts.mod error_handler.mod io_gcnmp.mod nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/thermalization.f90
