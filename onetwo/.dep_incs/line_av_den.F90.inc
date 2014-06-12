line_av_den.o: bicubic_spline.mod common_constants.mod cubic_spline.mod grid_class.mod io_gcnmp.mod mpi_data.mod nrtype.mod plasma_properties.mod solcon_gcnmp.mod vector_class.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/line_av_den.F90
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/line_av_den.F90
