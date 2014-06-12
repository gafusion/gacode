cray208.o: bicube.mod cer.mod constnts.mod contour.mod etc.mod extra.mod flags.mod flxav.mod geom.mod gpsi.mod ifs.mod io.mod limiter.mod machin.mod mesh.mod metrics.mod mhdbcdtn.mod mhdcom.mod mhdgrid.mod mhdpar.mod neo2dp.mod numbrs.mod param.mod psig.mod replace_imsl.mod rhog.mod shapctr.mod solcon.mod soln2d.mod soln.mod tfact.mod cray208.f
	$(COMPILE) -I$(NETCDF_INC) -DONETWO cray208.f
