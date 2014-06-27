zen_state.o: bd_condtn.mod common_constants.mod flags.mod fusion.mod io.mod ions.mod nrtype.mod nub2.mod numbrs.mod param.mod solcon.mod soln.mod verbose.mod zen_state.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo zen_state.mod | tr A-Z a-z | cut -d. -f1` zen_state.f90
zen_state.mod: bd_condtn.mod common_constants.mod flags.mod fusion.mod io.mod ions.mod nrtype.mod nub2.mod numbrs.mod param.mod solcon.mod soln.mod verbose.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO zen_state.f90
