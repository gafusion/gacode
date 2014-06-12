set_eqdsk_statefile_values.o: constnts.mod cubic_spline.mod flxav.mod gpsi.mod grid_class.mod mhdcom.mod nrtype.mod plasma_properties.mod psig.mod vector_class.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/set_eqdsk_statefile_values.f90
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/set_eqdsk_statefile_values.f90
