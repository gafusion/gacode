shot_info.o: nrtype.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/shot_info.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo shot_info.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/shot_info.f90
shot_info.mod: nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/shot_info.f90
