plasma_prop.o: ions_gcnmp.mod nrtype.mod vector_class.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/plasma_prop.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo plasma_properties.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/plasma_prop.f90
plasma_properties.mod: ions_gcnmp.mod nrtype.mod vector_class.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/plasma_prop.f90
