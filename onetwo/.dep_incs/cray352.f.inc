cray352.o: contour.mod ename.mod ext_prog_info.mod io.mod mhdpar.mod param.mod replace_imsl.mod solcon.mod cray352.f
	$(COMPILE) -I$(NETCDF_INC) -DONETWO cray352.f
