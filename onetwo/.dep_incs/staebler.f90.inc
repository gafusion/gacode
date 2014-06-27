staebler.o: param.mod staebler.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo staebler.mod | tr A-Z a-z | cut -d. -f1` staebler.f90
staebler.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO staebler.f90
