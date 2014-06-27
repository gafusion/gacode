cray501.o: adaptive.mod geom.mod ifs.mod mesh.mod mhdpar.mod neo2d.mod numbrs.mod param.mod replace_imsl.mod solcon.mod soln2d.mod soln.mod verbose.mod cray501.f
	$(COMPILE) -I$(NETCDF_INC) -DONETWO cray501.f
