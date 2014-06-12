mcgo.o: param.mod mcgo.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo mcgo.mod | tr A-Z a-z | cut -d. -f1` mcgo.f90
mcgo.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO mcgo.f90
