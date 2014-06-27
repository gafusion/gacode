neucg.o: param.mod neucg.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo neucg.mod | tr A-Z a-z | cut -d. -f1` neucg.f90
neucg.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO neucg.f90
