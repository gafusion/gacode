neutral_data.o: nrtype.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/neutral_data.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo neutral_data.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/neutral_data.f90
neutral_data.mod: nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/neutral_data.f90
