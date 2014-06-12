cray700.o: aid_newton.mod com.mod copy.mod flags.mod ifs.mod machin.mod mesh.mod nonlin.mod numbrs.mod param.mod solcon.mod tcoef.mod terminate.mod tordlrot.mod cray700.f
	$(COMPILE) -I$(NETCDF_INC) -DONETWO cray700.f
