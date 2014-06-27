nub.o: nrtype.mod param.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/nub.F90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo nub.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/nub.F90
nub.mod: nrtype.mod param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/nub.F90
