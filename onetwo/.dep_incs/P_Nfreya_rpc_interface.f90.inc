P_Nfreya_rpc_interface.o: common_constants.mod error_handler.mod ext_prog_info.mod file_proc.mod gcnmp_input.mod io_gcnmp.mod iterdbmd_gcnmp.mod iterdbmd.mod nfreya_namelist.mod nrtype.mod p_nfreya_12_interface.mod solcon_gcnmp.mod solcon.mod string_util.mod transp.mod P_Nfreya_rpc_interface.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo p_nfreya_rpc_interface.mod | tr A-Z a-z | cut -d. -f1` P_Nfreya_rpc_interface.f90
p_nfreya_rpc_interface.mod: common_constants.mod error_handler.mod ext_prog_info.mod file_proc.mod gcnmp_input.mod io_gcnmp.mod iterdbmd_gcnmp.mod iterdbmd.mod nfreya_namelist.mod nrtype.mod p_nfreya_12_interface.mod solcon_gcnmp.mod solcon.mod string_util.mod transp.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO P_Nfreya_rpc_interface.f90
