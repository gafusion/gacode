neo2d.o: param.mod neo2d.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo neo2d.mod | tr A-Z a-z | cut -d. -f1` neo2d.f90
neo2d.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO neo2d.f90
