mhdgrid.o: mhdpar.mod mhdgrid.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo mhdgrid.mod | tr A-Z a-z | cut -d. -f1` mhdgrid.f90
mhdgrid.mod: mhdpar.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO mhdgrid.f90
