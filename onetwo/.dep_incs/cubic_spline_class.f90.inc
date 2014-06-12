cubic_spline_class.o: common_constants.mod nrtype.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/cubic_spline_class.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo cubic_spline.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/cubic_spline_class.f90
cubic_spline.mod: common_constants.mod nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/cubic_spline_class.f90
