iterdbmd_gcnmp.o: nrtype.mod solcon_gcnmp.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/iterdbmd_gcnmp.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo iterdbmd_gcnmp.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/iterdbmd_gcnmp.f90
iterdbmd_gcnmp.mod: nrtype.mod solcon_gcnmp.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/iterdbmd_gcnmp.f90
