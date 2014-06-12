m_lines.o: m_lines.f90
	$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo m_lines.mod | tr A-Z a-z | cut -d. -f1` m_lines.f90
m_lines.mod:
	$(COMPILE) -I$(NETCDF_INC) -DONETWO m_lines.f90
