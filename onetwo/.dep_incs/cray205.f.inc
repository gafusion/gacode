cray205.o: bd_condtn.mod cer.mod constnts.mod contour.mod etc.mod extra.mod flxav.mod geom.mod io.mod ions.mod iterdbmd.mod machin.mod mesh.mod mhdpar.mod neut.mod nub2.mod nub.mod numbrs.mod param.mod psig.mod rhog.mod solcon.mod soln.mod sourc.mod tcoef.mod tordlrot.mod yoka.mod cray205.f
	$(COMPILE) -I$(NETCDF_INC) -DONETWO cray205.f
