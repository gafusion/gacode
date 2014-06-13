brainfuse_mod.o: string_util.mod brainfuse_mod.F90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo brainfuse_mod.mod | tr A-Z a-z | cut -d. -f1` brainfuse_mod.F90
brainfuse_mod.mod: string_util.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO brainfuse_mod.F90
