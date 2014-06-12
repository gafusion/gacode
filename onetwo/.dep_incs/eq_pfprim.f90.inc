eq_pfprim.o: constnts.mod geom.mod io.mod machin.mod mhdcom.mod numbrs.mod param.mod rhog.mod soln.mod eq_pfprim.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo eq_pfprim.mod | tr A-Z a-z | cut -d. -f1` eq_pfprim.f90
eq_pfprim.mod: constnts.mod geom.mod io.mod machin.mod mhdcom.mod numbrs.mod param.mod rhog.mod soln.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO eq_pfprim.f90
