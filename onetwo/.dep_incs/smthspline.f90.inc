smthspline.o: nrtype.mod smthspline.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo smthspline.mod | tr A-Z a-z | cut -d. -f1` smthspline.f90
smthspline.mod: nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO smthspline.f90
