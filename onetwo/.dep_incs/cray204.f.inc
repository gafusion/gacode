cray204.o: bd_condtn.mod constnts.mod contour.mod etc.mod extra.mod flags.mod flxav.mod fusion.mod geom.mod io.mod ions.mod iterdbmd.mod machin.mod mesh.mod mhdpar.mod neo2d.mod neut.mod nub2.mod nub.mod numbrs.mod param.mod psig.mod rhog.mod solcon.mod soln.mod sourc.mod tcoef.mod tordlrot.mod yoka.mod cray204.f
	$(COMPILE) -I$(NETCDF_INC) -DONETWO cray204.f
