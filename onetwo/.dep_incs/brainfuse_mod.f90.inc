brainfuse_mod.o: brainfuse_mod.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo brainfuse_mod.mod | tr A-Z a-z | cut -d. -f1` brainfuse_mod.f90
brainfuse_mod.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO brainfuse_mod.f90
