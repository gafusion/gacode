iterdbmd.o: nrtype.mod iterdbmd.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo iterdbmd.mod | tr A-Z a-z | cut -d. -f1` iterdbmd.f90
iterdbmd.mod: nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO iterdbmd.f90
