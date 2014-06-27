neo2dp.o: mhdpar.mod neo2dp.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo neo2dp.mod | tr A-Z a-z | cut -d. -f1` neo2dp.f90
neo2dp.mod: mhdpar.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO neo2dp.f90
