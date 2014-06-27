param.o: param.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo param.mod | tr A-Z a-z | cut -d. -f1` param.f90
param.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO param.f90
