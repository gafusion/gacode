fusion.o: param.mod fusion.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo fusion.mod | tr A-Z a-z | cut -d. -f1` fusion.f90
fusion.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO fusion.f90
