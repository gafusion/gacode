gcnmp_input.o: nrtype.mod param.mod solcon_gcnmp.mod gcnmp_input.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo gcnmp_input.mod | tr A-Z a-z | cut -d. -f1` gcnmp_input.f90
gcnmp_input.mod: nrtype.mod param.mod solcon_gcnmp.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO gcnmp_input.f90
