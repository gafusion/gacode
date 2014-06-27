load_12.o: fusion.mod machin.mod mesh.mod nub.mod numbrs.mod param.mod sourc.mod string_util.mod transp.mod verbose.mod load_12.f90
	$(COMPILE) -I$(NETCDF_INC) -DONETWO load_12.f90
