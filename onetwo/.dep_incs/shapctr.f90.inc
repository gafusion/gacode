shapctr.o: mhdpar.mod shapctr.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo shapctr.mod | tr A-Z a-z | cut -d. -f1` shapctr.f90
shapctr.mod: mhdpar.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO shapctr.f90
