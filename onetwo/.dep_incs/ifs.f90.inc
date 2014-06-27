ifs.o: mhdpar.mod param.mod ifs.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo ifs.mod | tr A-Z a-z | cut -d. -f1` ifs.f90
ifs.mod: mhdpar.mod param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO ifs.f90
