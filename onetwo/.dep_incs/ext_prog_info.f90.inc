ext_prog_info.o: io.mod mhdpar.mod ext_prog_info.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo ext_prog_info.mod | tr A-Z a-z | cut -d. -f1` ext_prog_info.f90
ext_prog_info.mod: io.mod mhdpar.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO ext_prog_info.f90
