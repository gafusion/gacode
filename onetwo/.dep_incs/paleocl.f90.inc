paleocl.o: aid_newton.mod constnts.mod extra.mod geom.mod io.mod ions.mod machin.mod mesh.mod numbrs.mod param.mod soln.mod sourc.mod tcoef.mod paleocl.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo paleocl.mod | tr A-Z a-z | cut -d. -f1` paleocl.f90
paleocl.mod: aid_newton.mod constnts.mod extra.mod geom.mod io.mod ions.mod machin.mod mesh.mod numbrs.mod param.mod soln.mod sourc.mod tcoef.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO paleocl.f90
