flx.o: param.mod flx.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo flx.mod | tr A-Z a-z | cut -d. -f1` flx.f90
flx.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO flx.f90
