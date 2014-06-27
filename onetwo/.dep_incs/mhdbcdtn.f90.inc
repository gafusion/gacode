mhdbcdtn.o: mhdpar.mod nrtype.mod mhdbcdtn.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo mhdbcdtn.mod | tr A-Z a-z | cut -d. -f1` mhdbcdtn.f90
mhdbcdtn.mod: mhdpar.mod nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO mhdbcdtn.f90
