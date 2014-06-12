nbnamelist.o: nrtype.mod numbrs.mod param.mod transp.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/nbnamelist.F90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo nbnamelist.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/nbnamelist.F90
nbnamelist.mod: nrtype.mod numbrs.mod param.mod transp.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/nbnamelist.F90
