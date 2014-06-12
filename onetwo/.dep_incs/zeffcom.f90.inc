zeffcom.o: nrtype.mod zeffcom.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo zeffcom.mod | tr A-Z a-z | cut -d. -f1` zeffcom.f90
zeffcom.mod: nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO zeffcom.f90
