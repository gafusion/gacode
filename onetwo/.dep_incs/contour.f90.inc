contour.o: contour.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo contour.mod | tr A-Z a-z | cut -d. -f1` contour.f90
contour.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO contour.f90
