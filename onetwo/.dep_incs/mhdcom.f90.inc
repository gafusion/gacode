mhdcom.o: mhdpar.mod mhdcom.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo mhdcom.mod | tr A-Z a-z | cut -d. -f1` mhdcom.f90
mhdcom.mod: mhdpar.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO mhdcom.f90
