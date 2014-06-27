cray311.o: geom.mod io.mod ions.mod mesh.mod neucg.mod neut.mod numbrs.mod param.mod replace_imsl.mod solcon.mod soln.mod sourc.mod verbose.mod cray311.f
	$(COMPILE) -I$(NETCDF_INC) -DONETWO cray311.f
