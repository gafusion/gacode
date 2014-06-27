modmmm7_1.o: nrtype.mod w20mod.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/modmmm7_1.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo modmmm7_1.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/modmmm7_1.f90
modmmm7_1.mod: nrtype.mod w20mod.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/modmmm7_1.f90
