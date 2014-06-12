xpload_module.o: xpload_module.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo xpload_module.mod | tr A-Z a-z | cut -d. -f1` xpload_module.f90
xpload_module.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO xpload_module.f90
