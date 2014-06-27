ode_pkg.o: ode_pkg.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo ode_pkg.mod | tr A-Z a-z | cut -d. -f1` ode_pkg.f90
ode_pkg.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO ode_pkg.f90
