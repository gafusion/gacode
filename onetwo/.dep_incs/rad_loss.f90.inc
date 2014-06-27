rad_loss.o: nrtype.mod rad_loss.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo rad_loss.mod | tr A-Z a-z | cut -d. -f1` rad_loss.f90
rad_loss.mod: nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO rad_loss.f90
