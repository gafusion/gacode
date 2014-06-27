pelcom.o: nrtype.mod pelcom.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo pelcom.mod | tr A-Z a-z | cut -d. -f1` pelcom.f90
pelcom.mod: nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO pelcom.f90
