soln.o: param.mod soln.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo soln.mod | tr A-Z a-z | cut -d. -f1` soln.f90
soln.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO soln.f90
