terminate.o: terminate.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo terminate.mod | tr A-Z a-z | cut -d. -f1` terminate.f90
terminate.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO terminate.f90
