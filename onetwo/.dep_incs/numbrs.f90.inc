numbrs.o: param.mod numbrs.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo numbrs.mod | tr A-Z a-z | cut -d. -f1` numbrs.f90
numbrs.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO numbrs.f90
