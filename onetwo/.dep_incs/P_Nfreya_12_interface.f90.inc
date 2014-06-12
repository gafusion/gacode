P_Nfreya_12_interface.o: common_constants.mod events.mod fusion.mod io_gcnmp.mod ions.mod neutral_beams.mod nrtype.mod nub2.mod nub3.mod nub.mod param.mod solcon.mod tordlrot.mod transp.mod P_Nfreya_12_interface.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo p_nfreya_12_interface.mod | tr A-Z a-z | cut -d. -f1` P_Nfreya_12_interface.f90
p_nfreya_12_interface.mod: common_constants.mod events.mod fusion.mod io_gcnmp.mod ions.mod neutral_beams.mod nrtype.mod nub2.mod nub3.mod nub.mod param.mod solcon.mod tordlrot.mod transp.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO P_Nfreya_12_interface.f90
