rhog.o: nrtype.mod param.mod rhog.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo rhog.mod | tr A-Z a-z | cut -d. -f1` rhog.f90
rhog.mod: nrtype.mod param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO rhog.f90
