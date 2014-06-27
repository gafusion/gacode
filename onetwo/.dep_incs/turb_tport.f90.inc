turb_tport.o: nrtype.mod solcon_gcnmp.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/turb_tport.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo turb_tport.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/turb_tport.f90
turb_tport.mod: nrtype.mod solcon_gcnmp.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/turb_tport.f90
