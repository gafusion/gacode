nubeam_driver.o: constnts.mod ename.mod ext_prog_info.mod geom.mod io.mod ions.mod machin.mod nbnamelist.mod neut.mod nub2.mod nub4.mod nub.mod numbrs.mod param.mod restore_12.mod solcon.mod soln.mod sourc.mod tordlrot.mod transp.mod nubeam_driver.f90
	$(COMPILE) -I$(NETCDF_INC) -DONETWO nubeam_driver.f90
