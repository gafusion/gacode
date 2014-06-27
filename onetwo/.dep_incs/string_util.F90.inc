string_util.o: error_handler.mod grid_class.mod io_gcnmp.mod mesh.mod nrtype.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/string_util.F90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo string_util.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/string_util.F90
string_util.mod: error_handler.mod grid_class.mod io_gcnmp.mod mesh.mod nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/string_util.F90
