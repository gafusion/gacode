cray203.o: bicube.mod constnts.mod gpsi.mod io.mod mhdcom.mod mhdgrid.mod mhdpar.mod nub.mod param.mod replace_imsl.mod cray203.f
	$(COMPILE) -I$(NETCDF_INC) -DONETWO cray203.f
