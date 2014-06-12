cray303.o: cer.mod fusion.mod geom.mod ions.mod machin.mod mesh.mod metrics.mod mhdpar.mod neo2d.mod nub2.mod numbrs.mod param.mod replace_imsl.mod soln.mod tfact.mod tordlrot.mod cray303.f
	$(COMPILE) -I$(NETCDF_INC) -DONETWO cray303.f
