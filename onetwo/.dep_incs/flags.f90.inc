flags.o: param.mod flags.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo flags.mod | tr A-Z a-z | cut -d. -f1` flags.f90
flags.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO flags.f90
