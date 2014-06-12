fdyimp.o: common_constants.mod grid_class.mod io_gcnmp.mod nrtype.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/fdyimp.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo fdyimp.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/fdyimp.f90
fdyimp.mod: common_constants.mod grid_class.mod io_gcnmp.mod nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/fdyimp.f90
