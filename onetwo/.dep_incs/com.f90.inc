com.o: com.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo com.mod | tr A-Z a-z | cut -d. -f1` com.f90
com.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO com.f90
