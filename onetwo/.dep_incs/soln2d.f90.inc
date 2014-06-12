soln2d.o: mhdpar.mod soln2d.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo soln2d.mod | tr A-Z a-z | cut -d. -f1` soln2d.f90
soln2d.mod: mhdpar.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO soln2d.f90
