metrics.o: mhdpar.mod metrics.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo metrics.mod | tr A-Z a-z | cut -d. -f1` metrics.f90
metrics.mod: mhdpar.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO metrics.f90
