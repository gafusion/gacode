vector_class.o: nrtype.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/vector_class.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo vector_class.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/vector_class.f90
vector_class.mod: nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/vector_class.f90
