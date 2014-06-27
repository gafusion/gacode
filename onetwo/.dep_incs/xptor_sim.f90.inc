xptor_sim.o: xptor_sim.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo xptor_sim.mod | tr A-Z a-z | cut -d. -f1` xptor_sim.f90
xptor_sim.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO xptor_sim.f90
