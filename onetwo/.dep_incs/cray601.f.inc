cray601.o: aid_newton.mod allocate_err.mod com.mod copy.mod flags.mod ifs.mod io.mod ions.mod machin.mod mesh.mod mhdpar.mod nonlin.mod numbrs.mod param.mod solcon.mod soln2d.mod soln.mod staebler.mod tcoef.mod terminate.mod tordlrot.mod verbose.mod cray601.f
	$(COMPILE) -I$(NETCDF_INC) -DONETWO cray601.f
