set_initial_profiles.o: bd_condtn.mod common_constants.mod flags.mod fusion.mod gpsi.mod io.mod ions.mod iterdbmd.mod mesh.mod nrtype.mod nub2.mod numbrs.mod param.mod psig.mod rf.mod solcon.mod soln.mod sourc.mod string_util.mod tordlrot.mod set_initial_profiles.f90
	$(COMPILE) -I$(NETCDF_INC) -DONETWO set_initial_profiles.f90
