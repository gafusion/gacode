solcon_gcnmp.o: common_constants.mod nrtype.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/solcon_gcnmp.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo solcon_gcnmp.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/solcon_gcnmp.f90
solcon_gcnmp.mod: common_constants.mod nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/solcon_gcnmp.f90
