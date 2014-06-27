fast_ion_data_gcnmp.o: nrtype.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/fast_ion_data_gcnmp.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo fast_ion_data_gcnmp.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/fast_ion_data_gcnmp.f90
fast_ion_data_gcnmp.mod: nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/fast_ion_data_gcnmp.f90
