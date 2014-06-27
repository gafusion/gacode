transp.o: beam_structure.mod fast_ion_diffusion.mod nrtype.mod transp.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo transp.mod | tr A-Z a-z | cut -d. -f1` transp.f90
transp.mod: beam_structure.mod fast_ion_diffusion.mod nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO transp.f90
