allocate_plasma_freq.o: grid_class.mod ions_gcnmp.mod nrtype.mod plasma_properties.mod vector_class.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/allocate_plasma_freq.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo pl_freq.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/allocate_plasma_freq.f90
pl_freq.mod: grid_class.mod ions_gcnmp.mod nrtype.mod plasma_properties.mod vector_class.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/allocate_plasma_freq.f90
