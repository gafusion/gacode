weiland.o: nrtype.mod param.mod weiland.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo weiland.mod | tr A-Z a-z | cut -d. -f1` weiland.f90
weiland.mod: nrtype.mod param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO weiland.f90
