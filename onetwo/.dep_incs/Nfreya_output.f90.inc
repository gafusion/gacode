Nfreya_output.o: common_constants.mod geom.mod io.mod mesh.mod nrtype.mod nub2.mod nub.mod numbrs.mod param.mod plasma_properties.mod p_nfreya_12_interface.mod transp.mod yoka.mod Nfreya_output.f90
	$(COMPILE) -I$(NETCDF_INC) -DONETWO Nfreya_output.f90
