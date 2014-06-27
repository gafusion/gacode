ename.o: solcon.mod yoka.mod ename.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo ename.mod | tr A-Z a-z | cut -d. -f1` ename.f90
ename.mod: solcon.mod yoka.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO ename.f90
