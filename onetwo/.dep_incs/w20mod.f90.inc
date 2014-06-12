w20mod.o: nrtype.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/w20mod.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo w20mod.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/w20mod.f90
w20mod.mod: nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/w20mod.f90
