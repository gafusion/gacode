gcnmp_version.o: nrtype.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/gcnmp_version.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo gcnmp_version.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/gcnmp_version.f90
gcnmp_version.mod: nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/gcnmp_version.f90
