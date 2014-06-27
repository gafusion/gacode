verbose.o: verbose.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo verbose.mod | tr A-Z a-z | cut -d. -f1` verbose.f90
verbose.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO verbose.f90
