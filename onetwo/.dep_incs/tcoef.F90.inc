tcoef.o: nf_param.mod nrtype.mod param.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/tcoef.F90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo tcoef.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/tcoef.F90
tcoef.mod: nf_param.mod nrtype.mod param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/tcoef.F90
