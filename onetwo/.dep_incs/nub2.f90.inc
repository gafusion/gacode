nub2.o: param.mod nub2.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo nub2.mod | tr A-Z a-z | cut -d. -f1` nub2.f90
nub2.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO nub2.f90
