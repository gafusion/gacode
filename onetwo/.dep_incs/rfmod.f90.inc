rfmod.o: ename.mod fusion.mod io.mod ions.mod nub2.mod nub.mod numbrs.mod param.mod solcon.mod soln.mod rfmod.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo rf.mod | tr A-Z a-z | cut -d. -f1` rfmod.f90
rf.mod: ename.mod fusion.mod io.mod ions.mod nub2.mod nub.mod numbrs.mod param.mod solcon.mod soln.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO rfmod.f90
