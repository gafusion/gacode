bicube.o: mhdpar.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/bicube.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo bicube.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/bicube.f90
bicube.mod: mhdpar.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/bicube.f90
