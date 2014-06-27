aid_newton.o: aid_newton.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo aid_newton.mod | tr A-Z a-z | cut -d. -f1` aid_newton.f90
aid_newton.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO aid_newton.f90
