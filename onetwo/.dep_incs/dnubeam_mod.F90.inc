dnubeam_mod.o: io_gcnmp.mod io.mod nrtype.mod param.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/dnubeam_mod.F90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo dnubeam_mod.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/dnubeam_mod.F90
dnubeam_mod.mod: io_gcnmp.mod io.mod nrtype.mod param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/dnubeam_mod.F90
