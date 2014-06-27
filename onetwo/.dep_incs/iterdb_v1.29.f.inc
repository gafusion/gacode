iterdb_v1.29.o: constnts.mod contour.mod etc.mod extra.mod flxav.mod geom.mod io.mod ions.mod iterdbmd.mod machin.mod mesh.mod mhdpar.mod neut.mod nub2.mod nub.mod numbrs.mod param.mod psig.mod rhog.mod solcon.mod soln.mod sourc.mod tcoef.mod tordlrot.mod yoka.mod iterdb_v1.29.f
	$(COMPILE) -I$(NETCDF_INC) -DONETWO iterdb_v1.29.f
