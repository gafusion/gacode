prl_mod.o: spec_kind_mod.mod prl_mod.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo prl_mod.mod | tr A-Z a-z | cut -d. -f1` prl_mod.f90
prl_mod.mod: spec_kind_mod.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO prl_mod.f90
