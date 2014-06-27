ions.o: param.mod ions.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo ions.mod | tr A-Z a-z | cut -d. -f1` ions.f90
ions.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO ions.f90
