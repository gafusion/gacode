neutral_beams.o: common_constants.mod ions_gcnmp.mod nf_param.mod nrtype.mod plasma_properties.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/neutral_beams.F90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo neutral_beams.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/neutral_beams.F90
neutral_beams.mod: common_constants.mod ions_gcnmp.mod nf_param.mod nrtype.mod plasma_properties.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/neutral_beams.F90
