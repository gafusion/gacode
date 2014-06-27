beam_structure.o: allocate_err.mod nrtype.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/beam_structure.F90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo beam_structure.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/beam_structure.F90
beam_structure.mod: allocate_err.mod nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/beam_structure.F90
