dnubeam_driver.o: beam_structure.mod constnts.mod dnubeam_mod.mod ename.mod ext_prog_info.mod extra.mod fusion.mod geom.mod io_gcnmp.mod io.mod ions.mod machin.mod nbnamelist.mod neut.mod nub2.mod nub4.mod nub.mod numbrs.mod param.mod plasma_properties.mod solcon.mod soln.mod sourc.mod tordlrot.mod transp.mod dnubeam_driver.f90
	$(COMPILE) -I$(NETCDF_INC) -DONETWO dnubeam_driver.f90
