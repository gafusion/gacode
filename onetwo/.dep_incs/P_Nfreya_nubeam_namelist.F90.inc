P_Nfreya_nubeam_namelist.o: dnubeam_mod.mod nrtype.mod numbrs.mod param.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/P_Nfreya_nubeam_namelist.F90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo p_nfreya_nubeam_namelist.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/P_Nfreya_nubeam_namelist.F90
p_nfreya_nubeam_namelist.mod: dnubeam_mod.mod nrtype.mod numbrs.mod param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/P_Nfreya_nubeam_namelist.F90
