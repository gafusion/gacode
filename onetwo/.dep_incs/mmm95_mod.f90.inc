mmm95_mod.o: mmm95_mod.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo mmm95_mod.mod | tr A-Z a-z | cut -d. -f1` mmm95_mod.f90
mmm95_mod.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO mmm95_mod.f90
