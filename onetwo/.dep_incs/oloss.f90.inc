oloss.o: param.mod oloss.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo oloss_dat.mod | tr A-Z a-z | cut -d. -f1` oloss.f90
oloss_dat.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO oloss.f90
