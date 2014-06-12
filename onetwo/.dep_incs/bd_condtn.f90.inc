bd_condtn.o: common_constants.mod nrtype.mod param.mod bd_condtn.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo bd_condtn.mod | tr A-Z a-z | cut -d. -f1` bd_condtn.f90
bd_condtn.mod: common_constants.mod nrtype.mod param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO bd_condtn.f90
