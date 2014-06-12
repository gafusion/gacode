kin_efit.o: bd_condtn.mod constnts.mod extra.mod fusion.mod geom.mod ions.mod machin.mod mesh.mod mhdcom.mod nrtype.mod nub2.mod nub.mod numbrs.mod param.mod psig.mod rhog.mod solcon.mod soln.mod sourc.mod yoka.mod kin_efit.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo kinetic_efit.mod | tr A-Z a-z | cut -d. -f1` kin_efit.f90
kinetic_efit.mod: bd_condtn.mod constnts.mod extra.mod fusion.mod geom.mod ions.mod machin.mod mesh.mod mhdcom.mod nrtype.mod nub2.mod nub.mod numbrs.mod param.mod psig.mod rhog.mod solcon.mod soln.mod sourc.mod yoka.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO kin_efit.f90
