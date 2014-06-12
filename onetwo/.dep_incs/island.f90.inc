island.o: nrtype.mod island.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo island.mod | tr A-Z a-z | cut -d. -f1` island.f90
island.mod: nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO island.f90
