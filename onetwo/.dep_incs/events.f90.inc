events.o: file_proc.mod io.mod mhdcom.mod mhdpar.mod nrtype.mod param.mod rf.mod solcon.mod tdem.mod events.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo events.mod | tr A-Z a-z | cut -d. -f1` events.f90
events.mod: file_proc.mod io.mod mhdcom.mod mhdpar.mod nrtype.mod param.mod rf.mod solcon.mod tdem.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO events.f90
