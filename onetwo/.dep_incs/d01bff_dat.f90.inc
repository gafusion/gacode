d01bff_dat.o: d01bff_dat.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo d01bff_dat.mod | tr A-Z a-z | cut -d. -f1` d01bff_dat.f90
d01bff_dat.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO d01bff_dat.f90
