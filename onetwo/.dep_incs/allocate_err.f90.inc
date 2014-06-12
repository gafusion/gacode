allocate_err.o: allocate_err.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo allocate_err.mod | tr A-Z a-z | cut -d. -f1` allocate_err.f90
allocate_err.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO allocate_err.f90
