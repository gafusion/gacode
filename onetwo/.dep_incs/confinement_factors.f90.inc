confinement_factors.o: bd_condtn.mod extra.mod fusion.mod ions.mod machin.mod nrtype.mod numbrs.mod soln.mod yoka.mod confinement_factors.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo cfactrs.mod | tr A-Z a-z | cut -d. -f1` confinement_factors.f90
cfactrs.mod: bd_condtn.mod extra.mod fusion.mod ions.mod machin.mod nrtype.mod numbrs.mod soln.mod yoka.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO confinement_factors.f90
