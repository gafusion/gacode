m_refsor.o: nrtype.mod /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/m_refsor.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo m_refsor.mod | tr A-Z a-z | cut -d. -f1` /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/m_refsor.f90
m_refsor.mod: nrtype.mod
	$(COMPILE) -I$(NETCDF_INC) -DONETWO /u/smithsp/source/gacode_onetwo/onetwo/shared_modules/m_refsor.f90
