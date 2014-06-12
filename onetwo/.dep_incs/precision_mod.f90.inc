precision_mod.o: precision_mod.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo precision_mod.mod | tr A-Z a-z | cut -d. -f1` precision_mod.f90
precision_mod.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO precision_mod.f90
