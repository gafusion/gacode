etc.o: nrtype.mod param.mod etc.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo etc.mod | tr A-Z a-z | cut -d. -f1` etc.f90
etc.mod: nrtype.mod param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO etc.f90
