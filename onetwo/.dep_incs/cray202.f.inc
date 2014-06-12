cray202.o: bicube.mod constnts.mod contour.mod gpsi.mod io.mod limiter.mod mhdbcdtn.mod mhdcom.mod mhdgrid.mod mhdpar.mod numbrs.mod param.mod replace_imsl.mod rhog.mod toq_12_interface.mod cray202.f
	$(COMPILE) -I$(NETCDF_INC) -DONETWO cray202.f
