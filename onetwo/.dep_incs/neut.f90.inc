neut.o: param.mod neut.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo neut.mod | tr A-Z a-z | cut -d. -f1` neut.f90
neut.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO neut.f90
