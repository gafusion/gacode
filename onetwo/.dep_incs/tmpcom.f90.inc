tmpcom.o: param.mod tmpcom.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo tmpcom.mod | tr A-Z a-z | cut -d. -f1` tmpcom.f90
tmpcom.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO tmpcom.f90
