nubeam_12.o: transp.mod nubeam_12.f
	$(COMPILE) -I$(NETCDF_INC) -DONETWO nubeam_12.f
