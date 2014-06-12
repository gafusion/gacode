cray404.o: flx.mod geom.mod ions.mod machin.mod mesh.mod monitr.mod neo2d.mod numbrs.mod param.mod soln.mod weiland.mod cray404.f
	$(COMPILE) -I$(NETCDF_INC) -DONETWO cray404.f
