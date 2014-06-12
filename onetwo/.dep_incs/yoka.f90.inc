yoka.o: param.mod yoka.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo yoka.mod | tr A-Z a-z | cut -d. -f1` yoka.f90
yoka.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO yoka.f90
