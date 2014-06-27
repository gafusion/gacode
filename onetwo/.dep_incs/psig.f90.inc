psig.o: mhdpar.mod nrtype.mod psig.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo psig.mod | tr A-Z a-z | cut -d. -f1` psig.f90
psig.mod: mhdpar.mod nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO psig.f90
