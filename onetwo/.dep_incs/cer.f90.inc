cer.o: mhdpar.mod nrtype.mod param.mod cer.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo cer.mod | tr A-Z a-z | cut -d. -f1` cer.f90
cer.mod: mhdpar.mod nrtype.mod param.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO cer.f90
