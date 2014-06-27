nub3.o: param.mod nub3.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo nub3.mod | tr A-Z a-z | cut -d. -f1` nub3.f90
nub3.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO nub3.f90
