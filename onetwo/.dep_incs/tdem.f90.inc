tdem.o: contour.mod limiter.mod mhdpar.mod param.mod tdem.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo tdem.mod | tr A-Z a-z | cut -d. -f1` tdem.f90
tdem.mod: contour.mod limiter.mod mhdpar.mod param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO tdem.f90
