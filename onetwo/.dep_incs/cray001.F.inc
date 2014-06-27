cray001.o: adaptive.mod constnts.mod ext_prog_info.mod io.mod mcgo.mod param.mod staebler.mod cray001.F
	$(COMPILE) -I$(NETCDF_INC) -DONETWO cray001.F
