machin.o: param.mod machin.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo machin.mod | tr A-Z a-z | cut -d. -f1` machin.f90
machin.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO machin.f90
