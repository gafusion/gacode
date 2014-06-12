replace_imsl.o: bicube.mod bicubic_spline.mod cubic_spline.mod error_handler.mod io_gcnmp.mod m_refsor.mod nrtype.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/replace_imsl.F90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo replace_imsl.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/replace_imsl.F90
replace_imsl.mod: bicube.mod bicubic_spline.mod cubic_spline.mod error_handler.mod io_gcnmp.mod m_refsor.mod nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/replace_imsl.F90
