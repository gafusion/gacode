cray304.o: cer.mod extra.mod fusion.mod geom.mod ions.mod machin.mod mesh.mod metrics.mod mhdpar.mod neo2d.mod nub2.mod numbrs.mod param.mod soln.mod tcoef.mod tordlrot.mod cray304.f
	$(COMPILE) -I$(NETCDF_INC) -DONETWO cray304.f
