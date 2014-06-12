nf_param.o: /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/nf_param.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo nf_param.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/nf_param.f90
nf_param.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/nf_param.f90
