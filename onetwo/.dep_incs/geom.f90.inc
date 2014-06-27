geom.o: param.mod geom.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo geom.mod | tr A-Z a-z | cut -d. -f1` geom.f90
geom.mod: param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO geom.f90
