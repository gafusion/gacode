glf23.o: param.mod glf23.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo glf23.mod | tr A-Z a-z | cut -d. -f1` glf23.f90
glf23.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO glf23.f90
