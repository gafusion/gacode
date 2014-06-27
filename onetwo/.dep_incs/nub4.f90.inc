nub4.o: param.mod nub4.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo nub4.mod | tr A-Z a-z | cut -d. -f1` nub4.f90
nub4.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO nub4.f90
