cray403.o: numbrs.mod param.mod soln.mod cray403.f
	$(COMPILE) -I$(NETCDF_INC) -DONETWO cray403.f
