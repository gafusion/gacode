restore_12.o: nrtype.mod param.mod restore_12.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo restore_12.mod | tr A-Z a-z | cut -d. -f1` restore_12.f90
restore_12.mod: nrtype.mod param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO restore_12.f90
