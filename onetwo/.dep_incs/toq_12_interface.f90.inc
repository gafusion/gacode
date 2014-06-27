toq_12_interface.o: ext_prog_info.mod io.mod mhdpar.mod toq_12_interface.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo toq_12_interface.mod | tr A-Z a-z | cut -d. -f1` toq_12_interface.f90
toq_12_interface.mod: ext_prog_info.mod io.mod mhdpar.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO toq_12_interface.f90
