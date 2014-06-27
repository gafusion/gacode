periodic_table_mod.o: periodic_table_mod.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo periodic_table_mod.mod | tr A-Z a-z | cut -d. -f1` periodic_table_mod.f90
periodic_table_mod.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO periodic_table_mod.f90
