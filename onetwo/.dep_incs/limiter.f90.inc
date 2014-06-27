limiter.o: limiter.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo limiter.mod | tr A-Z a-z | cut -d. -f1` limiter.f90
limiter.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO limiter.f90
