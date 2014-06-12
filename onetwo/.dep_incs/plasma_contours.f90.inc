plasma_contours.o: bicubic_spline.mod common_constants.mod error_handler.mod grid_class.mod io_gcnmp.mod mpi_data.mod nrtype.mod plasma_properties.mod vector_class.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/plasma_contours.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo plasma_contours.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/plasma_contours.f90
plasma_contours.mod: bicubic_spline.mod common_constants.mod error_handler.mod grid_class.mod io_gcnmp.mod mpi_data.mod nrtype.mod plasma_properties.mod vector_class.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/plasma_contours.f90
