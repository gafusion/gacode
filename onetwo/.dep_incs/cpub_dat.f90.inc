cpub_dat.o: cpub_dat.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo cpub_dat.mod | tr A-Z a-z | cut -d. -f1` cpub_dat.f90
cpub_dat.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO cpub_dat.f90
