tglf_input_defaults.o: common_constants.mod nrtype.mod solcon_gcnmp.mod tglfin.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/tglf_input_defaults.f90
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/tglf_input_defaults.f90
