mmm_input_defaults.o: common_constants.mod error_handler.mod io_gcnmp.mod mmm_in.mod modmmm7_1.mod nrtype.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/mmm_input_defaults.F90
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/mmm_input_defaults.F90
