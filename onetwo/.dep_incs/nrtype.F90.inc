nrtype.o: /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/nrtype.F90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo nrtype.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/nrtype.F90
nrtype.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/nrtype.F90
