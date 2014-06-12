tspline.o: error_handler.mod io_gcnmp.mod nrtype.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/tspline.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo tension_spline.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/tspline.f90
tension_spline.mod: error_handler.mod io_gcnmp.mod nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/tspline.f90
