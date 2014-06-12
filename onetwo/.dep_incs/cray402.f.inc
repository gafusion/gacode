cray402.o: colrate.mod common_constants.mod constnts.mod contour.mod flags.mod fusion.mod geom.mod ions.mod limiter.mod machin.mod mesh.mod mhdcom.mod mhdpar.mod neut.mod nub2.mod nub.mod numbrs.mod param.mod solcon.mod soln.mod sourc.mod tdem.mod tmpcom.mod tordlrot.mod verbose.mod cray402.f
	$(COMPILE) -I$(NETCDF_INC) -DONETWO cray402.f
