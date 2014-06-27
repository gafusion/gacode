tfact.o: param.mod tfact.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo tfact.mod | tr A-Z a-z | cut -d. -f1` tfact.f90
tfact.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO tfact.f90
