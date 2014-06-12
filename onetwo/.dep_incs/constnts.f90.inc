constnts.o: constnts.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo constnts.mod | tr A-Z a-z | cut -d. -f1` constnts.f90
constnts.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO constnts.f90
