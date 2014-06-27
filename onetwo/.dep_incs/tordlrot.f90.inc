tordlrot.o: param.mod tordlrot.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo tordlrot.mod | tr A-Z a-z | cut -d. -f1` tordlrot.f90
tordlrot.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO tordlrot.f90
