nubeam_sub.o: nubeam_sub.f90
	$(COMPILE) -I$(NETCDF_INC) -DONETWO nubeam_sub.f90
