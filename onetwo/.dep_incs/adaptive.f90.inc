adaptive.o: param.mod adaptive.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo adaptive.mod | tr A-Z a-z | cut -d. -f1` adaptive.f90
adaptive.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO adaptive.f90
