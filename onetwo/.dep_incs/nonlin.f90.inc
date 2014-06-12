nonlin.o: param.mod nonlin.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo nonlin.mod | tr A-Z a-z | cut -d. -f1` nonlin.f90
nonlin.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO nonlin.f90
