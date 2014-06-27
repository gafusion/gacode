io.o: param.mod io.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo io.mod | tr A-Z a-z | cut -d. -f1` io.f90
io.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO io.f90
