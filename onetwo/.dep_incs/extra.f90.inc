extra.o: param.mod extra.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo extra.mod | tr A-Z a-z | cut -d. -f1` extra.f90
extra.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO extra.f90
