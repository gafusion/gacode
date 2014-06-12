monitr.o: param.mod monitr.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo monitr.mod | tr A-Z a-z | cut -d. -f1` monitr.f90
monitr.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO monitr.f90
