check_en_bc.o: bd_condtn.mod common_constants.mod io.mod mesh.mod nrtype.mod numbrs.mod param.mod check_en_bc.f90
	$(COMPILE) -I$(NETCDF_INC) -DONETWO check_en_bc.f90
