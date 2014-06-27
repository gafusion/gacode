load_restart_12.o: geom.mod numbrs.mod restore_12.mod soln.mod sourc.mod tordlrot.mod load_restart_12.f90
	$(COMPILE) -I$(NETCDF_INC) -DONETWO load_restart_12.f90
