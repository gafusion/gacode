forcebal_data_mod.o: spec_kind_mod.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/forcebal_data_mod.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo forcebal_data_mod.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/forcebal_data_mod.f90
forcebal_data_mod.mod: spec_kind_mod.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/forcebal_data_mod.f90
