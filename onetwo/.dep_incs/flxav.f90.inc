flxav.o: flxav.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo flxav.mod | tr A-Z a-z | cut -d. -f1` flxav.f90
flxav.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO flxav.f90
