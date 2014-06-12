mixcom.o: param.mod mixcom.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo mixcom.mod | tr A-Z a-z | cut -d. -f1` mixcom.f90
mixcom.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO mixcom.f90
