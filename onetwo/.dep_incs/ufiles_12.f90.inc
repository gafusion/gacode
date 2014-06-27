ufiles_12.o: beam_structure.mod fast_ion_diffusion.mod fusion.mod geom.mod io.mod ions.mod nbnamelist.mod nub3.mod nub.mod numbrs.mod param.mod periodic_table_mod.mod p_nfreya_12_interface.mod solcon.mod string_util.mod tordlrot.mod transp.mod ufiles_12.f90
	$(COMPILE) -I$(NETCDF_INC) -DONETWO ufiles_12.f90
