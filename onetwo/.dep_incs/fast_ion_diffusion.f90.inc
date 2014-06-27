fast_ion_diffusion.o: fast_ion_diffusion.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo fast_ion_diffusion.mod | tr A-Z a-z | cut -d. -f1` fast_ion_diffusion.f90
fast_ion_diffusion.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO fast_ion_diffusion.f90
