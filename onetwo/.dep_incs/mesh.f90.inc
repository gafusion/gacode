mesh.o: param.mod mesh.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo mesh.mod | tr A-Z a-z | cut -d. -f1` mesh.f90
mesh.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO mesh.f90
