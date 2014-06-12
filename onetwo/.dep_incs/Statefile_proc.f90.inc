Statefile_proc.o: file_proc.mod gcnmp_input.mod iterdbmd_gcnmp.mod iterdbmd.mod mhdcom.mod nrtype.mod numbrs.mod set_12_gcnmp_vars.mod solcon.mod soln2d.mod Statefile_proc.f90
	$(COMPILE) -I$(NETCDF_INC) -DONETWO Statefile_proc.f90
