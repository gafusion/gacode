copy.o: copy.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo copy.mod | tr A-Z a-z | cut -d. -f1` copy.f90
copy.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO copy.f90
