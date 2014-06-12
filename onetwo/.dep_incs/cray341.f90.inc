cray341.o: events.mod fusion.mod geom.mod io.mod ions.mod machin.mod mesh.mod neut.mod nub2.mod nub.mod numbrs.mod param.mod pelcom.mod prl_mod.mod solcon.mod soln.mod sourc.mod spec_kind_mod.mod cray341.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo pellet_mod.mod | tr A-Z a-z | cut -d. -f1` cray341.f90
pellet_mod.mod: events.mod fusion.mod geom.mod io.mod ions.mod machin.mod mesh.mod neut.mod nub2.mod nub.mod numbrs.mod param.mod pelcom.mod prl_mod.mod solcon.mod soln.mod sourc.mod spec_kind_mod.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO cray341.f90
