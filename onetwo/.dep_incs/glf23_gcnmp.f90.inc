glf23_gcnmp.o: nrtype.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/glf23_gcnmp.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo glf23_gcnmp.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/glf23_gcnmp.f90
glf23_gcnmp.mod: nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/glf23_gcnmp.f90
