nbi_restart.o: ename.mod fusion.mod mhdcom.mod nub2.mod nub.mod numbrs.mod param.mod solcon.mod soln.mod sourc.mod tordlrot.mod transp.mod verbose.mod nbi_restart.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo nbi_restart.mod | tr A-Z a-z | cut -d. -f1` nbi_restart.f90
nbi_restart.mod: ename.mod fusion.mod mhdcom.mod nub2.mod nub.mod numbrs.mod param.mod solcon.mod soln.mod sourc.mod tordlrot.mod transp.mod verbose.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO nbi_restart.f90
