echdat_module.o: error_handler.mod io_gcnmp.mod io.mod nrtype.mod param.mod rf.mod echdat_module.F90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo echdat_module.mod | tr A-Z a-z | cut -d. -f1` echdat_module.F90
echdat_module.mod: error_handler.mod io_gcnmp.mod io.mod nrtype.mod param.mod rf.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO echdat_module.F90
