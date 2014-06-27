file_proc.o: echdat_module.mod io.mod nrtype.mod solcon.mod file_proc.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo file_proc.mod | tr A-Z a-z | cut -d. -f1` file_proc.f90
file_proc.mod: echdat_module.mod io.mod nrtype.mod solcon.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO file_proc.f90
