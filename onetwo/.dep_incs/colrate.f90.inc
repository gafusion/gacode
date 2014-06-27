colrate.o: colrate.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo colrate.mod | tr A-Z a-z | cut -d. -f1` colrate.f90
colrate.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO colrate.f90
