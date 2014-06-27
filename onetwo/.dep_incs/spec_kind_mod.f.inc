spec_kind_mod.o: spec_kind_mod.f
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo spec_kind_mod.mod | tr A-Z a-z | cut -d. -f1` spec_kind_mod.f
spec_kind_mod.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO spec_kind_mod.f
