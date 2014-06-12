bicubic_spline.o: error_handler.mod io_gcnmp.mod nrtype.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/bicubic_spline.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo bicubic_spline.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/bicubic_spline.f90
bicubic_spline.mod: error_handler.mod io_gcnmp.mod nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/bicubic_spline.f90
