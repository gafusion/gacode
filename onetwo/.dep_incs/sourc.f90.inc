sourc.o: param.mod sourc.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo sourc.mod | tr A-Z a-z | cut -d. -f1` sourc.f90
sourc.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO sourc.f90
