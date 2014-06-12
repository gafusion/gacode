gpsi.o: nrtype.mod param.mod gpsi.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo gpsi.mod | tr A-Z a-z | cut -d. -f1` gpsi.f90
gpsi.mod: nrtype.mod param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO gpsi.f90
