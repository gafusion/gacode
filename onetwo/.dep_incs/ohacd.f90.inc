ohacd.o: bd_condtn.mod constnts.mod extra.mod fusion.mod geom.mod io.mod ions.mod machin.mod mesh.mod mixcom.mod nub2.mod numbrs.mod param.mod rf.mod solcon.mod soln2d.mod soln.mod sourc.mod tcoef.mod tfact.mod tordlrot.mod verbose.mod ohacd.f90
	$(COMPILE) -I$(NETCDF_INC) -DONETWO ohacd.f90
