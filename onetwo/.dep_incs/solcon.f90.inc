solcon.o: param.mod solcon.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo solcon.mod | tr A-Z a-z | cut -d. -f1` solcon.f90
solcon.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO solcon.f90
