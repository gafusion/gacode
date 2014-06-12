mhdpar.o: param.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/mhdpar.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo mhdpar.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/mhdpar.f90
mhdpar.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/mhdpar.f90
